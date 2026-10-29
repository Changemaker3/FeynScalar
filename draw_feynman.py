#!/usr/bin/env python3
"""Draw Feynman diagrams exported by FeynScalar.wl.

    python3 draw_feynman.py diagrams.json -o diagrams.pdf
    python3 draw_feynman.py diagrams.json -o diagrams.png --cols 3 --momenta --amplitude
    python3 draw_feynman.py diagrams.json -o all.pdf --separate "figs/diagram_{id}.pdf"

The JSON file is written by ExportDiagrams[...] (or DrawDiagrams[...]) in FeynScalar.wl.
Needs numpy and matplotlib; scipy is used for the layout when it is installed.

Line styles: real scalars are solid when the theory has only real fields, otherwise dashed;
complex scalars are solid with an arrow showing the particle flow (--real-style overrides).
"""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import defaultdict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.mathtext import MathTextParser  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

try:
    from scipy.optimize import minimize as _scipy_minimize
except Exception:  # scipy is optional
    _scipy_minimize = None

FIELD_COLORS = ["#1a1a1a", "#1f5fbf", "#c0392b", "#1e8449", "#7d3c98", "#b9770e"]
VERTEX_COLORS = ["#1a1a1a", "#e67e22", "#8e44ad", "#16a085", "#c0392b", "#7f8c8d"]
DASHES = (0, (4.5, 2.8))
LINE_WIDTH = 1.6
_PARSER = MathTextParser("path")


# --------------------------------------------------------------------------- TeX helpers

def clean_tex(s: str) -> str:
    """Turn Mathematica TeXForm output into something matplotlib's mathtext accepts."""
    for a, b in [("\\text{", "\\mathrm{"), ("\\left(", "("), ("\\right)", ")"),
                 ("\\left[", "["), ("\\right]", "]"), ("\\left\\{", "\\{"),
                 ("\\right\\}", "\\}"), ("\\left.", ""), ("\\right.", ""), ("{}^", "^")]:
        s = s.replace(a, b)
    return s.strip()


def as_math(s: str, fallback: str | None = None) -> str | None:
    """'$s$' if mathtext can render it, else the fallback."""
    t = "$" + clean_tex(s) + "$"
    try:
        _PARSER.parse(t)
        return t
    except Exception:
        return fallback


# --------------------------------------------------------------------------- geometry

def bezier(a, b, h, n=80):
    """Cubic Bezier from a to b whose midpoint sits h (signed) to the left of a->b.
    h = +-|b-a|/2 gives a semicircle."""
    d = b - a
    L = np.linalg.norm(d)
    nrm = np.array([-d[1], d[0]]) / max(L, 1e-12)
    alpha = 0.0 if abs(h) > 0.45 * L else 1.0 / 3.0   # semicircle, or a flatter lens
    c1 = a + alpha * d + (4.0 / 3.0) * h * nrm
    c2 = b - alpha * d + (4.0 / 3.0) * h * nrm
    tt = np.linspace(0.0, 1.0, n)[:, None]
    return ((1 - tt) ** 3) * a + 3 * ((1 - tt) ** 2) * tt * c1 + 3 * (1 - tt) * tt ** 2 * c2 + tt ** 3 * b


def _orient(p, q, r):
    return (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])


def polyline_crossings(A, B) -> int:
    """Number of proper crossings between two polylines (arrays of points)."""
    if len(A) < 2 or len(B) < 2:
        return 0
    P1, P2 = A[:-1][:, None, :], A[1:][:, None, :]
    Q1, Q2 = B[:-1][None, :, :], B[1:][None, :, :]
    d1, d2 = _orient(Q1, Q2, P1), _orient(Q1, Q2, P2)
    d3, d4 = _orient(P1, P2, Q1), _orient(P1, P2, Q2)
    return int(np.count_nonzero((d1 * d2 < -1e-12) & (d3 * d4 < -1e-12)))


def point_polyline_distance(p, A) -> float:
    a, b = A[:-1], A[1:]
    ab = b - a
    t = np.clip(np.einsum("ij,ij->i", p - a, ab) / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-12), 0, 1)
    return float(np.min(np.linalg.norm(p - (a + t[:, None] * ab), axis=1)))


def best_gaps(angles, k):
    """k directions that bisect the largest angular gaps between the given angles."""
    if not angles:
        return [math.pi / 2 + 2 * math.pi * j / k for j in range(k)]
    out = []
    angs = sorted(a % (2 * math.pi) for a in angles)
    for _ in range(k):
        gaps = []
        for i, a in enumerate(angs):
            b = angs[(i + 1) % len(angs)] + (2 * math.pi if i == len(angs) - 1 else 0.0)
            gaps.append((b - a, a))
        width = max(g[0] for g in gaps)
        cands = [a + w / 2 for w, a in gaps if w > width - 1e-3]
        best = max(cands, key=lambda c: math.sin(c))  # prefer pointing up
        out.append(best)
        angs = sorted(angs + [best % (2 * math.pi)])
    return out


def arc_offsets(m, mode="round"):
    """Bulges of m parallel lines, in units of the vertex distance."""
    if m == 1:
        return [0.0]
    if mode == "lens":
        return [0.22, -0.22] if m == 2 else list(np.linspace(0.32, -0.32, m))
    if m <= 3:
        return list(np.linspace(0.5, -0.5, m))
    return list(np.linspace(0.55, -0.55, m))


# --------------------------------------------------------------------------- layout

class Layout:
    """Positions of the internal vertices. Incoming legs end on the left, outgoing legs on the
    right; internal vertices minimise a spring + repulsion energy. Many random starts are tried
    and ranked by the crossings and near-misses of the lines as they are actually drawn."""

    def __init__(self, diag, restarts=48, seed=0):
        self.V = len(diag["vertices"])
        self.ext = diag["externals"]
        self.edges = diag["edges"]
        self.E = len(self.ext)
        self.rng = np.random.default_rng(seed)

        self.groups = defaultdict(list)      # (a, b), a < b  ->  edge indices
        self.loop_edges = defaultdict(list)  # vertex -> self-loop edge indices
        for idx, e in enumerate(self.edges):
            a, b = e["tail"] - 1, e["head"] - 1
            if a == b:
                self.loop_edges[a].append(idx)
            else:
                self.groups[(min(a, b), max(a, b))].append(idx)

        d_io, d_side = self._distances()
        self.W = max(3.0, 0.95 * d_io)
        self.L0 = self.W / max(d_io, 2)
        self.H = max(1.0, 0.6 * self.L0 * (d_side - 2))
        self.anchor = self._anchors()
        self.mode = "round"
        self.pos = self._optimise(restarts)

    # external endpoints -------------------------------------------------
    def _anchors(self):
        ins = [e for e in self.ext if e["direction"] == "in"]
        outs = [e for e in self.ext if e["direction"] == "out"]
        pos = {}
        if ins and outs:
            for group, x in ((ins, 0.0), (outs, self.W)):
                n = len(group)
                h = self.H * (1.0 + 0.3 * max(0, n - 2))
                ys = [0.0] if n == 1 else np.linspace(h, -h, n)
                for e, y in zip(group, ys):
                    pos[e["id"]] = np.array([x, float(y)])
        else:  # only incoming or only outgoing legs: put them on a circle
            group = ins or outs
            n = len(group)
            R = 0.5 * self.W
            for k, e in enumerate(group):
                ang = math.pi - 2 * math.pi * k / n
                pos[e["id"]] = np.array([0.5 * self.W + R * math.cos(ang), R * math.sin(ang)])
        return pos

    def _distances(self):
        """Longest graph distance between an incoming and an outgoing leg, and between two legs
        on the same side."""
        V = self.V
        adj = defaultdict(set)
        for (a, b) in self.groups:
            adj[a].add(b)
            adj[b].add(a)
        for k, e in enumerate(self.ext):
            x, v = V + k, e["vertex"] - 1
            adj[x].add(v)
            adj[v].add(x)

        def bfs(src):
            dist, frontier = {src: 0}, [src]
            while frontier:
                nxt = []
                for x in frontier:
                    for y in adj[x]:
                        if y not in dist:
                            dist[y] = dist[x] + 1
                            nxt.append(y)
                frontier = nxt
            return dist

        ins = [V + k for k, e in enumerate(self.ext) if e["direction"] == "in"]
        outs = [V + k for k, e in enumerate(self.ext) if e["direction"] == "out"]
        if not ins or not outs:
            return 2 + max(0, V - 1), 2
        d_io, d_side = 2, 2
        for src in ins + outs:
            dist = bfs(src)
            same = ins if src in ins else outs
            other = outs if src in ins else ins
            d_io = max([d_io] + [dist[x] for x in other if x in dist])
            d_side = max([d_side] + [dist[x] for x in same if x in dist and x != src])
        return d_io, d_side

    # curves as drawn ---------------------------------------------------
    def curves(self, P, n=80, mode=None):
        """Every line of the diagram as a polyline: dicts with kind, index, pts, nodes."""
        V, out = self.V, []
        mode = mode or self.mode
        for k, e in enumerate(self.ext):
            v = e["vertex"] - 1
            out.append({"kind": "ext", "index": k, "pts": np.linspace(self.anchor[e["id"]], P[v], max(n // 4, 8)),
                        "nodes": {V + k, v}, "offset": 0.0, "pair": None})
        for (a, b), idxs in self.groups.items():
            A, B = P[a], P[b]
            d = np.linalg.norm(B - A)
            for off, idx in zip(arc_offsets(len(idxs), mode), idxs):
                out.append({"kind": "edge", "index": idx, "pts": bezier(A, B, off * d, n),
                            "nodes": {a, b}, "offset": off, "pair": (a, b)})
        r = 0.32 * self.L0
        for v, idxs in self.loop_edges.items():
            angles = []
            for c in out:
                if v in c["nodes"]:
                    pts = c["pts"] if np.linalg.norm(c["pts"][0] - P[v]) < 1e-9 else c["pts"][::-1]
                    q = pts[max(1, len(pts) // 5)] - P[v]
                    angles.append(math.atan2(q[1], q[0]))
            for ang, idx in zip(best_gaps(angles, len(idxs)), idxs):
                cen = P[v] + r * np.array([math.cos(ang), math.sin(ang)])
                th = np.linspace(ang + math.pi, ang + 3 * math.pi, max(n, 24))
                out.append({"kind": "loop", "index": idx, "nodes": {v}, "offset": 0.0, "pair": None,
                            "angle": ang, "center": cen, "radius": r,
                            "pts": np.column_stack([cen[0] + r * np.cos(th), cen[1] + r * np.sin(th)])})
        return out

    # energy ---------------------------------------------------------------
    def _setup(self):
        V, E, L0 = self.V, self.E, self.L0
        self.F = np.array([self.anchor[e["id"]] for e in self.ext], dtype=float).reshape(E, 2)
        si, sj, sl, sk = [], [], [], []
        for (a, b), idxs in self.groups.items():
            m = len(idxs)
            si.append(a)
            sj.append(b)
            sl.append(L0 * (1.0 if m == 1 else 1.0 + 0.12 * m))
            sk.append(float(m))
        for k, e in enumerate(self.ext):
            si.append(e["vertex"] - 1)
            sj.append(V + k)
            sl.append(0.85 * L0)
            sk.append(1.0)
        self.SI, self.SJ = np.array(si, int), np.array(sj, int)
        self.SL, self.SK = np.array(sl, float), np.array(sk, float)
        ra, rb, rw = [], [], []
        nl = {v: len(x) for v, x in self.loop_edges.items()}
        for i in range(V):
            for j in range(i + 1, V + E):
                w = 0.35 * L0 ** 2 * (1.0 + 0.6 * nl.get(i, 0) + (0.6 * nl.get(j, 0) if j < V else 0.0))
                ra.append(i)
                rb.append(j)
                rw.append(w)
        self.RA, self.RB, self.RW = np.array(ra, int), np.array(rb, int), np.array(rw, float)

    def _energy(self, x):
        V = self.V
        P = np.vstack([x.reshape(V, 2), self.F])
        g = np.zeros_like(P)
        e = 0.0
        if len(self.SI):
            d = P[self.SI] - P[self.SJ]
            r = np.linalg.norm(d, axis=1) + 1e-9
            e += np.sum(self.SK * (r - self.SL) ** 2)
            c = (2 * self.SK * (r - self.SL) / r)[:, None] * d
            np.add.at(g, self.SI, c)
            np.add.at(g, self.SJ, -c)
        if len(self.RA):
            d = P[self.RA] - P[self.RB]
            r = np.linalg.norm(d, axis=1) + 1e-9
            e += np.sum(self.RW / r)
            c = (-self.RW / r ** 3)[:, None] * d
            np.add.at(g, self.RA, c)
            np.add.at(g, self.RB, -c)
        # keep vertices between the incoming and the outgoing legs
        X = P[:V]
        lo, hi, ymax = 0.08 * self.W, 0.92 * self.W, self.H + 0.9 * self.L0
        for k, (low, high) in ((0, (lo, hi)), (1, (-ymax, ymax))):
            below = np.minimum(X[:, k] - low, 0.0)
            above = np.maximum(X[:, k] - high, 0.0)
            e += 8.0 * np.sum(below ** 2 + above ** 2)
            g[:V, k] += 16.0 * (below + above)
        return e, g[:V].ravel()

    def _minimise(self, x0):
        if _scipy_minimize is not None:
            res = _scipy_minimize(self._energy, x0, jac=True, method="L-BFGS-B", options={"maxiter": 400})
            return res.x, float(res.fun)
        x, step = x0.copy(), 0.05
        e, g = self._energy(x)
        for _ in range(1500):
            xn = x - step * g
            en, gn = self._energy(xn)
            if en < e:
                x, e, g, step = xn, en, gn, step * 1.2
            else:
                step *= 0.5
                if step < 1e-8:
                    break
        return x, e

    def _score(self, x, energy):
        """(score, bulge mode): round bubbles unless flatter lenses avoid crossings."""
        modes = ["round", "lens"] if any(len(v) > 1 for v in self.groups.values()) else ["round"]
        return min((self._score_mode(x, energy, md) + (0.5 if md == "lens" else 0.0), md) for md in modes)

    def _score_mode(self, x, energy, mode):
        V, L0 = self.V, self.L0
        P = x.reshape(V, 2)
        allpos = {i: P[i] for i in range(V)}
        allpos.update({V + k: self.anchor[e["id"]] for k, e in enumerate(self.ext)})
        cs = self.curves(P, n=16, mode=mode)
        cross = 0
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                A, B = cs[i], cs[j]
                if A["pair"] is not None and A["pair"] == B["pair"]:
                    continue
                pa, pb = A["pts"], B["pts"]
                for node in A["nodes"] & B["nodes"]:   # ignore the shared endpoint
                    c = allpos[node]
                    pa = pa[np.linalg.norm(pa - c, axis=1) > 0.12 * L0]
                    pb = pb[np.linalg.norm(pb - c, axis=1) > 0.12 * L0]
                cross += polyline_crossings(pa, pb)
        near = 0
        for v in range(V):
            for c in cs:
                if v not in c["nodes"] and point_polyline_distance(P[v], c["pts"]) < 0.2 * L0:
                    near += 1
            for w in range(v + 1, V):
                if np.linalg.norm(P[v] - P[w]) < 0.4 * L0:
                    near += 1
        return 10.0 * cross + 4.0 * near + energy

    def _optimise(self, restarts):
        V = self.V
        if V == 0:
            return np.zeros((0, 2))
        self._setup()
        starts = []
        guess = np.zeros((V, 2))
        for v in range(V):
            pts = [self.anchor[e["id"]] for e in self.ext if e["vertex"] - 1 == v]
            guess[v] = (0.6 * np.mean(pts, axis=0) + np.array([0.2 * self.W, 0.0])) if pts else \
                np.array([0.5 * self.W, 0.0])
        starts.append(guess + 0.05 * self.rng.standard_normal((V, 2)))
        for _ in range(restarts):
            starts.append(np.column_stack([
                self.rng.uniform(0.15 * self.W, 0.85 * self.W, V),
                self.rng.uniform(-0.9 * self.H, 0.9 * self.H, V)]))
        best, best_score = None, math.inf
        for x0 in starts:
            x, e = self._minimise(x0.ravel())
            sc, md = self._score(x, e)
            if sc < best_score - 1e-9:
                best, best_score, self.mode = x, sc, md
        return best.reshape(V, 2)


# --------------------------------------------------------------------------- styles

def field_styles(fields, real_style="auto"):
    only_real = all(f["type"] == "real" for f in fields)
    styles = {}
    for k, f in enumerate(fields):
        color = FIELD_COLORS[k % len(FIELD_COLORS)] if len(fields) > 1 else FIELD_COLORS[0]
        if f["type"] == "complex":
            ls = "solid"
        elif real_style == "solid" or (real_style == "auto" and only_real):
            ls = "solid"
        else:
            ls = DASHES
        styles[f["index"]] = {"color": color, "ls": ls, "arrow": f["type"] == "complex"}
    return styles


def draw_arrow(ax, pt, tangent, color, size):
    u = tangent / max(np.linalg.norm(tangent), 1e-12)
    n = np.array([-u[1], u[0]])
    tip, base = pt + 0.55 * size * u, pt - 0.45 * size * u
    ax.add_patch(Polygon([tip, base + 0.42 * size * n, base - 0.42 * size * n], closed=True,
                         facecolor=color, edgecolor=color, lw=0.5, zorder=5))


def draw_curve(ax, pts, style, arrow_forward=None, size=0.13):
    ax.plot(pts[:, 0], pts[:, 1], color=style["color"], ls=style["ls"], lw=LINE_WIDTH,
            solid_capstyle="round", zorder=2)
    if style["arrow"] and arrow_forward is not None:
        i = len(pts) // 2
        tangent = pts[min(i + 1, len(pts) - 1)] - pts[max(i - 1, 0)]
        draw_arrow(ax, pts[i], tangent if arrow_forward else -tangent, style["color"], size)


# --------------------------------------------------------------------------- one diagram

def draw_diagram(ax, diag, data, styles, opts, show_process=False):
    lay = Layout(diag, seed=diag["id"])
    V, L0, P = lay.V, lay.L0, lay.pos
    fields = {f["index"]: f for f in data["fields"]}
    multi_terms = len(data["terms"]) > 1
    show_particles = opts.particles if opts.particles is not None else \
        (len(fields) > 1 or any(ff["type"] == "complex" for ff in fields.values()))
    arrow_size = 0.14 * L0
    centroid = P.mean(axis=0) if V else np.zeros(2)
    xs, ys = [], []

    def track(pts):
        pts = np.atleast_2d(pts)
        xs.extend(pts[:, 0])
        ys.extend(pts[:, 1])

    for c in lay.curves(P):
        pts = c["pts"]
        track(pts)
        if c["kind"] == "ext":
            e = diag["externals"][c["index"]]
            fwd = {"to_vertex": True, "from_vertex": False}.get(e["arrow"])
            draw_curve(ax, pts, styles[e["field"]], fwd, arrow_size)
            a = pts[0]
            mom = as_math(e["momentum_tex"], "")
            label = mom
            if show_particles:
                name = clean_tex(fields[e["field"]]["tex"])
                if e["anti"]:
                    name = "\\bar{" + name + "}"
                label = as_math(name + "(" + e["momentum_tex"] + ")", mom)
            if a[0] <= 1e-9:
                ax.text(a[0] - 0.08, a[1], label, ha="right", va="center", fontsize=10)
                track([a[0] - 0.6, a[1]])
            elif a[0] >= lay.W - 1e-9:
                ax.text(a[0] + 0.08, a[1], label, ha="left", va="center", fontsize=10)
                track([a[0] + 0.6, a[1]])
            else:
                d = a - np.array([0.5 * lay.W, 0.0])
                d = d / max(np.linalg.norm(d), 1e-9)
                ax.text(*(a + 0.2 * d), label, ha="center", va="center", fontsize=10)
                track(a + 0.45 * d)
            continue
        e = diag["edges"][c["index"]]
        if c["kind"] == "edge":
            fwd = (e["tail"] - 1 == c["pair"][0]) if e["directed"] else None
        else:
            fwd = True if e["directed"] else None
        draw_curve(ax, pts, styles[e["field"]], fwd, arrow_size)
        if opts.momenta:
            lab = as_math(e["momentum_tex"])
            if not lab:
                continue
            if c["kind"] == "loop":
                tip = c["center"] + 1.45 * c["radius"] * np.array([math.cos(c["angle"]), math.sin(c["angle"])])
                ax.text(*tip, lab, fontsize=7, ha="center", va="center", color="#555555")
                track(tip)
            else:
                i = len(pts) // 2
                mid, tang = pts[i], pts[i + 1] - pts[i - 1]
                nrm = np.array([-tang[1], tang[0]]) / max(np.linalg.norm(tang), 1e-9)
                if abs(c["offset"]) < 1e-9:
                    side = 1.0 if np.dot(nrm, mid - centroid) >= 0 else -1.0
                else:
                    side = 1.0 if c["offset"] > 0 else -1.0
                ax.text(*(mid + side * 0.17 * L0 * nrm), lab, fontsize=7, ha="center", va="center",
                        color="#555555")

    for k, vv in enumerate(diag["vertices"]):
        col = VERTEX_COLORS[(vv["term"] - 1) % len(VERTEX_COLORS)] if multi_terms else VERTEX_COLORS[0]
        ax.plot(*P[k], "o", ms=5.5, color=col, zorder=6)
        track(P[k])

    S = diag["symmetry_factor"]
    caption = f"({diag['id']})  {diag['label']}" + ("" if S == 1 else f"    S = {S}")
    lines = [caption]
    if show_process and diag.get("process_tex"):
        proc = as_math(diag["process_tex"])
        if proc:
            lines.insert(0, proc)
    if opts.amplitude:
        amp = as_math(diag["amplitude_tex"])
        if amp and len(diag["amplitude_tex"]) < 160:
            lines.append(amp)
    pad = 0.35
    x0, x1 = min(xs) - pad, max(xs) + pad
    y0, y1 = min(ys) - pad, max(ys) + pad
    y = y0
    for k, txt in enumerate(lines):
        ax.text(0.5 * (x0 + x1), y, txt, ha="center", va="top", fontsize=9.5 if k < 2 else 8.5)
        y -= 0.34
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0 - 0.1 - 0.34 * len(lines), y1)
    ax.set_aspect("equal")
    ax.axis("off")


# --------------------------------------------------------------------------- figure

def legend_handles(data, styles, diagrams):
    used_fields = set()
    used_terms = set()
    for d in diagrams:
        used_fields.update(e["field"] for e in d["edges"])
        used_fields.update(e["field"] for e in d["externals"])
        used_terms.update(v["term"] for v in d["vertices"])
    handles = []
    for f in data["fields"]:
        if f["index"] not in used_fields:
            continue
        st = styles[f["index"]]
        name = as_math(f["tex"], f["name"])
        mass = as_math(f["mass_tex"], f["mass_tex"])
        extra = ", arrow = particle flow" if f["type"] == "complex" else ""
        handles.append(Line2D([0], [0], color=st["color"], ls=st["ls"], lw=LINE_WIDTH,
                              label=f"{name}  ({f['type']}, mass {mass}{extra})"))
    if len(data["terms"]) > 1:
        for tm in data["terms"]:
            if tm["index"] not in used_terms:
                continue
            col = VERTEX_COLORS[(tm["index"] - 1) % len(VERTEX_COLORS)]
            label = "vertex " + (as_math(tm["factor_tex"], tm["factor_tex"]) or "")
            term = as_math(tm.get("term_tex", ""))
            if term:
                label += "   from " + term
            handles.append(Line2D([0], [0], marker="o", ls="none", color=col, ms=6, label=label))
    return handles


def render(data, diagrams, out_paths, opts, title=True):
    styles = field_styles(data["fields"], opts.real_style)
    n = len(diagrams)
    cols = opts.cols or min(4, n)
    rows = math.ceil(n / cols)
    processes = {d.get("process_tex", data.get("title", "")) for d in diagrams}
    several = len(processes) > 1      # diagrams of several processes: label each panel
    fig, axes = plt.subplots(rows, cols, figsize=(3.4 * cols, 2.9 * rows + 0.6), squeeze=False)
    for ax, d in zip(axes.flat, diagrams):
        draw_diagram(ax, d, data, styles, opts, show_process=several)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    handles = legend_handles(data, styles, diagrams)
    title = title and not several and bool(processes - {""})
    if title:
        t = next(iter(processes))
        fig.suptitle(as_math(t, t), fontsize=13, y=0.995)
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 3), frameon=False,
                   fontsize=9, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.06 if handles else 0, 1, 0.96 if title else 1))
    for path in out_paths:
        d = os.path.dirname(os.path.abspath(path))
        os.makedirs(d, exist_ok=True)
        fig.savefig(path, dpi=opts.dpi, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json", help="diagram file written by ExportDiagrams[] in FeynScalar.wl")
    ap.add_argument("-o", "--output", action="append", default=None,
                    help="output file (.pdf, .png, .svg); may be given several times")
    ap.add_argument("--cols", type=int, default=None, help="diagrams per row (default: up to 4)")
    ap.add_argument("--momenta", action="store_true", help="label internal lines with their momenta")
    ap.add_argument("--amplitude", action="store_true", help="print i M under each diagram")
    ap.add_argument("--no-particles", dest="particles", action="store_false", default=None,
                    help="label external legs with momenta only")
    ap.add_argument("--real-style", choices=["auto", "solid", "dashed"], default="auto",
                    help="line style of real scalars")
    ap.add_argument("--separate", metavar="PATTERN", default=None,
                    help='also save each diagram alone, e.g. "figs/diagram_{id}.pdf"')
    ap.add_argument("--only", type=str, default=None, help='draw only these diagram ids, e.g. "1,3,5"')
    ap.add_argument("--dpi", type=int, default=200)
    opts = ap.parse_args(argv)

    with open(opts.json, encoding="utf-8") as fh:
        data = json.load(fh)
    diagrams = data["diagrams"]
    if opts.only:
        keep = {int(x) for x in opts.only.split(",") if x.strip()}
        diagrams = [d for d in diagrams if d["id"] in keep]
    if not diagrams:
        raise SystemExit("no diagrams to draw")
    outs = opts.output or [os.path.splitext(opts.json)[0] + ".pdf"]
    render(data, diagrams, outs, opts)
    if opts.separate:
        for d in diagrams:
            opts_one = argparse.Namespace(**vars(opts))
            opts_one.cols = 1
            render(data, [d], [opts.separate.replace("{id}", str(d["id"]))], opts_one, title=False)
    for path in outs:
        print(path)


if __name__ == "__main__":
    main()
