(* ::Package:: *)

(* ::Title:: *)
(*FeynScalar*)

(* ::Text:: *)
(*Feynman diagrams, symmetry factors and scattering amplitudes for theories of*)
(*real and complex scalar fields with polynomial interactions.*)
(**)
(*Conventions*)
(*  metric (+,-,-,-);  propagator  i/(p^2 - m^2 + i eps)*)
(*  L = L_free + Lint;  vertex factor = i * dLint/d(fields)   ( -lambda/4! phi^4  ->  -i lambda )*)
(*  i M = sum of connected, amputated diagrams,   <f|S-1|i> = (2 pi)^4 delta^4(P_f - P_i) i M*)
(*  momenta: p[1..nin] incoming, p[nin+1..n] outgoing (physical momenta, sum in = sum out)*)
(*  2 -> 2:  s = (p1+p2)^2,  t = (p1-p3)^2,  u = (p1-p4)^2*)
(*  one-loop integrals in LoopTools normalisation:*)
(*     Integral d^4l/(2pi)^4  1/[((l+r1)^2-m1^2)...((l+rn)^2-mn^2)]  =  i/(16 pi^2) X0,   X0 = A0, B0, C0, D0*)
(*  DeltaUV = 1/eps - EulerGamma + Log[4 pi]  with  D = 4 - 2 eps  (set DeltaUV -> 0 for MS-bar)*)
(*Diagrams are drawn by draw_feynman.py (Python + matplotlib), called from DrawDiagrams[].*)


BeginPackage["FeynScalar`"];

RealScalar::usage = "RealScalar[\[Phi], m] declares a real scalar field \[Phi] with mass m. RealScalar[\[Phi]] is massless.";
ComplexScalar::usage = "ComplexScalar[\[Chi], M] declares a complex scalar field \[Chi] with mass M. In the Lagrangian write the conjugate field as Conjugate[\[Chi]] (typed \[Chi] Esc co Esc) or use Abs[\[Chi]]^2. ComplexScalar[\[Chi], M, \[Chi]b] also accepts the symbol \[Chi]b for the conjugate field.";
Anti::usage = "Anti[\[Chi]] is the antiparticle of the complex field \[Chi] in a process. Conjugate[\[Chi]] and the conjugate alias (if declared) mean the same.";
ScalarTheory::usage = "ScalarTheory[{fields...}, Lint] defines a theory from field declarations and the interaction Lagrangian Lint (a polynomial in the fields).";
FeynmanRules::usage = "FeynmanRules[theory] shows the propagators and vertex factors of a theory.";
FeynmanDiagrams::usage = "FeynmanDiagrams[theory, {in...} -> {out...}, opts] generates all Feynman diagrams of a process with their symmetry factors and amplitudes.
Options:
  \"Loops\" -> L or {L1, L2, ...}      loop order (default: tree level)
  \"Order\" -> n or {n1, n2, ...}      number of vertices (order in perturbation theory)
  \"Order\" -> {g -> n, ...}           exact powers of the couplings (couplings not listed: power 0)
  \"Diagrams\" -> \"Amputated\" (default) | \"1PI\" | \"Connected\"
  \"Tadpoles\" -> True (default) | False
  \"OnShell\" -> Automatic | True | False   (Automatic: on shell, except for 2-point functions and \"Connected\")";
FeynmanDiagram::usage = "FeynmanDiagram[assoc] is one diagram. d[\"iM\"], d[\"M\"], d[\"SymmetryFactor\"], d[\"Integrand\"], d[\"Edges\"], d[\"Momenta\"], d[\"Label\"] ... read its data; d[\"Properties\"] lists them.";
DiagramTable::usage = "DiagramTable[diagrams] tabulates the diagrams with symmetry factors and i M. DiagramTable[diagrams, \"Integrand\"] shows the loop integrands instead.";
Amplitude::usage = "Amplitude[diagrams] gives the total amplitude M (where i M is the sum of the diagrams). Amplitude[diagrams, \"iM\"] gives i M.";
PVExplicit::usage = "PVExplicit[expr] writes A0 in closed form and B0, C0, D0 as Feynman-parameter integrals (use Activate to evaluate them).";
EvaluateLoops::usage = "EvaluateLoops[expr, {sym -> value, ...}] substitutes the values and evaluates A0, B0, C0, D0 numerically. Options \"DeltaUV\" -> value and \"Mu\" -> value fix the UV pole and the scale (default: kept as symbols; \"DeltaUV\" -> 0 is MS-bar).";
ExportDiagrams::usage = "ExportDiagrams[\"file.json\", diagrams] writes the diagrams as JSON for draw_feynman.py.";
DrawDiagrams::usage = "DrawDiagrams[diagrams, opts] draws the diagrams with Python/matplotlib and returns the picture.
Options: \"File\" -> \"out.pdf\" (also save; .pdf/.png/.svg), \"Columns\" -> n, \"ShowMomenta\" -> True, \"ShowAmplitude\" -> True, \"Separate\" -> \"dir/diagram_{id}.pdf\", \"LineStyle\" -> Automatic | \"Solid\" | \"Dashed\".";

A0::usage = "A0[m^2] is the one-loop scalar tadpole integral (LoopTools convention).";
B0::usage = "B0[p^2, m1^2, m2^2] is the one-loop scalar 2-point integral (LoopTools convention).";
C0::usage = "C0[p1^2, p2^2, (p1+p2)^2, m1^2, m2^2, m3^2] is the one-loop scalar 3-point integral (LoopTools convention).";
D0::usage = "D0[p1^2, p2^2, p3^2, p4^2, (p1+p2)^2, (p2+p3)^2, m1^2, m2^2, m3^2, m4^2] is the one-loop scalar 4-point integral (LoopTools convention).";
ScalarNPoint::usage = "ScalarNPoint[{r1, ..., rn}, {m1^2, ..., mn^2}] is the one-loop scalar n-point integral (n >= 5) with denominators (l + ri)^2 - mi^2, normalised like A0..D0.";
LoopIntegral::usage = "LoopIntegral[{l[1], l[2], ...}, {Den[q1, m1], ...}] is the multi-loop integral \[Integral] \[Product]d^4l/(2\[Pi])^4 \[Product] 1/(q^2 - m^2 + i\[Epsilon]).";
Den::usage = "Den[q, m] is the propagator denominator q^2 - m^2.";
SP::usage = "SP[a, b] is the Minkowski product a\[CenterDot]b.";
p::usage = "p[i] is the momentum of external particle i: 1..nin incoming, then outgoing.";
l::usage = "l[i] is the i-th loop momentum.";
s::usage = "Mandelstam s = (p1 + p2)^2.";
t::usage = "Mandelstam t = (p1 - p3)^2.";
u::usage = "Mandelstam u = (p1 - p4)^2.";
DeltaUV::usage = "DeltaUV = 1/\[Epsilon] - EulerGamma + Log[4\[Pi]] (D = 4 - 2\[Epsilon]), the UV pole of dimensional regularisation. Set it to 0 for MS-bar.";
MuR::usage = "MuR is the renormalisation scale \[Mu] of dimensional regularisation.";
$FeynScalarPython::usage = "$FeynScalarPython is the Python executable used by DrawDiagrams (found automatically if not set).";
$FeynScalarDirectory::usage = "$FeynScalarDirectory is the folder that contains FeynScalar.wl and draw_feynman.py.";

ScalarTheory::field = "`1` is not a field declaration. Use RealScalar[\[Phi], m] or ComplexScalar[\[Chi], M].";
ScalarTheory::dup = "A field is declared twice.";
ScalarTheory::poly = "The interaction Lagrangian `1` is not a polynomial in the declared fields.";
ScalarTheory::noint = "The interaction Lagrangian contains no interaction terms.";
FeynmanDiagrams::ext = "`1` is not a particle of this theory.";
FeynmanDiagrams::empty = "The process has no external particles.";
FeynmanDiagrams::order = "The number of vertices is not bounded (the theory has 1- or 2-point vertices). Give \"Order\" -> n, the number of vertices.";
FeynmanDiagrams::dtype = "\"Diagrams\" must be \"Amputated\", \"1PI\" or \"Connected\", not `1`.";
FeynmanDiagrams::none = "No diagrams contribute at this order.";
DrawDiagrams::nopy = "No Python with numpy and matplotlib was found. Set $FeynScalarPython = \"/path/to/python3\".";
DrawDiagrams::pyerr = "draw_feynman.py failed: `1`";

Begin["`Private`"];

$FeynScalarDirectory = DirectoryName[$InputFileName];


(* ::Section:: *)
(*Display*)


SetAttributes[SP, Orderless];
Format[p[i_Integer]] := Subscript[p, i];
Format[l[i_Integer]] := Subscript["\[ScriptL]", i];
Format[SP[a_, a_]] := Superscript[a, 2];
Format[SP[a_, b_]] := CenterDot[a, b];
Format[Den[q_, 0]] := HoldForm[q^2];
Format[Den[q_, m_]] := HoldForm[q^2 - m^2];
Format[LoopIntegral[ls_List, dens_List]] :=
  Row[{"\[Integral]", Row[Row[{Superscript["d", 4], #}] & /@ ls, " "], "/", Superscript["(2\[Pi])", 4 Length[ls]],
    " ", HoldForm[1/#] &[Times @@ dens]}];

Format[ScalarTheory[a_Association]] :=
  Row[{"ScalarTheory[", Row[fieldShow /@ a["Fields"], ", "], ";  Lint = ", a["Lint"], "]"}];
fieldShow[f_] := Row[{f["Symbol"], "(", If[f["Type"] == "Real", "real", "complex"], ", m = ", f["Mass"], ")"}];

Format[FeynmanDiagram[a_Association]] :=
  Row[{"FeynmanDiagram[#", a["Id"], ": ", a["Label"], ", S = ", a["SymmetryFactor"], "]"}];
FeynmanDiagram[a_Association][key_String] := Which[key === "Properties", Keys[a], True, a[key]];


(* ::Section:: *)
(*Theory*)


ScalarTheory[f : (_RealScalar | _ComplexScalar), lint_] := ScalarTheory[{f}, lint];
ScalarTheory[fields_List, lint_] := Catch[ScalarTheory[makeTheory[fields, lint]], $fsTag];

parseField[RealScalar[f : (_Symbol | _Subscript)], k_] := parseField[RealScalar[f, 0], k];
parseField[RealScalar[f : (_Symbol | _Subscript), m_], k_] :=
  <|"Symbol" -> f, "Type" -> "Real", "Mass" -> m, "Index" -> k, "Alias" -> None|>;
parseField[ComplexScalar[f : (_Symbol | _Subscript)], k_] := parseField[ComplexScalar[f, 0], k];
parseField[ComplexScalar[f : (_Symbol | _Subscript), m_], k_] :=
  <|"Symbol" -> f, "Type" -> "Complex", "Mass" -> m, "Index" -> k, "Alias" -> None|>;
parseField[ComplexScalar[f : (_Symbol | _Subscript), m_, a : (_Symbol | _Subscript)], k_] :=
  <|"Symbol" -> f, "Type" -> "Complex", "Mass" -> m, "Index" -> k, "Alias" -> a|>;
parseField[x_, _] := (Message[ScalarTheory::field, x]; Throw[$Failed, $fsTag]);

(* Internal variables: X[k, 0] real field k; X[k, 1] complex field k; X[k, -1] its conjugate *)
fieldRules[f_] := With[{k = f["Index"], x = f["Symbol"]},
   If[f["Type"] == "Real",
    {x -> X[k, 0]},
    Join[{Conjugate[x] -> X[k, -1], Abs[x] -> Sqrt[X[k, 1] X[k, -1]]},
     If[f["Alias"] =!= None, {f["Alias"] -> X[k, -1]}, {}], {x -> X[k, 1]}]]];

makeTheory[fields_List, lint_] :=
 Module[{flds, st, conj, vars, disp, poly, cr, terms},
  flds = MapIndexed[parseField[#1, First[#2]] &, fields];
  If[! DuplicateFreeQ[#Symbol & /@ flds], Message[ScalarTheory::dup]; Throw[$Failed, $fsTag]];
  (* slot types: one per real field, two (field, conjugate) per complex field *)
  st = Join @@ (If[#Type == "Real", {{#Index, 0}}, {{#Index, 1}, {#Index, -1}}] & /@ flds);
  conj = Position[st, {#[[1]], -#[[2]]}, {1}][[1, 1]] & /@ st;
  vars = X @@@ st;
  disp = Function[{k, sg}, With[{f = flds[[k]]},
       Which[sg >= 0, f["Symbol"], f["Alias"] =!= None, f["Alias"], True, Conjugate[f["Symbol"]]]]] @@@ st;
  poly = Expand[lint /. Join @@ (fieldRules /@ flds)];
  If[! PolynomialQ[poly, vars], Message[ScalarTheory::poly, lint]; Throw[$Failed, $fsTag]];
  cr = Select[CoefficientRules[poly, vars], Total[First[#]] > 0 &];
  If[cr === {}, Message[ScalarTheory::noint]; Throw[$Failed, $fsTag]];
  cr = SortBy[cr, {Total[First[#]] &, -First[#] &}];
  terms = MapIndexed[
    Function[{r, i}, <|"Index" -> First[i], "Exponents" -> r[[1]], "Coefficient" -> r[[2]],
      "Degree" -> Total[r[[1]]], "Factor" -> Simplify[I r[[2]] Times @@ Factorial[r[[1]]]],
      "Monomial" -> Times @@ (disp^r[[1]]),
      "Legs" -> Flatten[MapThread[ConstantArray, {disp, r[[1]]}]]|>], cr];
  <|"Fields" -> flds, "SlotTypes" -> st, "Conj" -> conj, "RealQ" -> (#Type == "Real" & /@ flds),
    "Terms" -> terms, "Lint" -> lint, "DisplayVars" -> disp|>]


FeynmanRules[ScalarTheory[th_Association]] :=
 Module[{props, verts},
  props = Function[f, {
      Row[{f["Symbol"], If[f["Type"] == "Complex", "   complex (arrow = particle flow)", "   real"]}],
      If[f["Mass"] === 0, HoldForm[I/("k"^2 + I "\[Epsilon]")],
       HoldForm[I/("k"^2 - mm^2 + I "\[Epsilon]")] /. mm -> f["Mass"]]}] /@ th["Fields"];
  verts = Function[tm, {tm["Coefficient"] tm["Monomial"], Row[tm["Legs"], " "], tm["Factor"]}] /@ th["Terms"];
  Column[{
    Style["Propagators", Bold],
    Grid[Prepend[props, Style[#, Italic] & /@ {"field", "propagator"}], Frame -> All, Alignment -> Left, Spacings -> {2, 1}],
    "",
    Style["Vertices   (vertex factor = i \[Times] coefficient \[Times] \[Product] n!;  momentum is conserved at each vertex)", Bold],
    Grid[Prepend[verts, Style[#, Italic] & /@ {"term in Lint", "legs", "vertex factor"}], Frame -> All, Alignment -> Left, Spacings -> {2, 1}]
    }, Spacings -> 1]]


(* ::Section:: *)
(*Process and kinematics*)


parseProcess[th_, in_ -> out_] :=
 Module[{ins = Flatten[{in}], outs = Flatten[{out}]},
  If[Length[ins] + Length[outs] == 0, Message[FeynmanDiagrams::empty]; Throw[$Failed, $fsTag]];
  Join[MapIndexed[parseExt[th, #1, "In", First[#2]] &, ins],
   MapIndexed[parseExt[th, #1, "Out", Length[ins] + First[#2]] &, outs]]]

parseExt[th_, x_, dir_, j_] :=
 Module[{flds = th["Fields"], anti = False, f = x, k, fd, sg},
  Which[
   MatchQ[x, _Anti | _Conjugate], f = First[x]; anti = True,
   MemberQ[DeleteCases[#Alias & /@ flds, None], x], f = SelectFirst[flds, #Alias === x &]["Symbol"]; anti = True];
  k = First[FirstPosition[#Symbol & /@ flds, f, {Missing[]}, {1}, Heads -> False]];
  If[MissingQ[k], Message[FeynmanDiagrams::ext, x]; Throw[$Failed, $fsTag]];
  fd = flds[[k]];
  If[fd["Type"] == "Real", anti = False];
  (* incoming particle / outgoing antiparticle contract with chi; the others with chi^* *)
  sg = If[fd["Type"] == "Real", 0, If[Xor[dir == "In", anti], 1, -1]];
  <|"Index" -> j, "Field" -> k, "Symbol" -> fd["Symbol"], "Anti" -> anti, "Direction" -> dir,
    "SlotType" -> {k, sg}, "Slot" -> Position[th["SlotTypes"], {k, sg}, {1}][[1, 1]],
    "Momentum" -> p[j], "Mass" -> fd["Mass"]|>]

extTeX[e_] := If[e["Anti"], "\\bar{" <> texStr[e["Symbol"]] <> "}", texStr[e["Symbol"]]];
processTeX[exts_] :=
  StringRiffle[extTeX /@ Select[exts, #Direction == "In" &], "\\,"] <> " \\to " <>
   StringRiffle[extTeX /@ Select[exts, #Direction == "Out" &], "\\,"];

(* Replacement rules for products of external momenta *)
kinRules[exts_, onShell_] :=
 Module[{n = Length[exts], m2, sg, sij},
  m2[i_] := If[onShell, exts[[i, "Mass"]]^2, SP[p[i], p[i]]];
  sg[i_] := If[exts[[i, "Direction"]] == "In", 1, -1];
  sij[i_, j_] := If[n == 3, m2[First[Complement[{1, 2, 3}, {i, j}]]],
    Switch[Sort[{i, j}], {1, 2} | {3, 4}, s, {1, 3} | {2, 4}, t, {1, 4} | {2, 3}, u]];
  <|"Rules" -> Which[
     n == 1, {},
     n == 2, If[onShell, {SP[p[1], p[1]] -> exts[[1, "Mass"]]^2}, {}],
     n == 3 || n == 4,
     Join[If[onShell, Table[SP[p[i], p[i]] -> m2[i], {i, n}], {}],
      Flatten[Table[SP[p[i], p[j]] -> sg[i] sg[j] (sij[i, j] - m2[i] - m2[j])/2, {i, n}, {j, i + 1, n}]]],
     True, If[onShell, Table[SP[p[i], p[i]] -> m2[i], {i, n}], {}]],
    (* s + t + u = sum of the external masses squared *)
    "Sum" -> If[n == 4, Sum[m2[i], {i, 4}], None]|>]

(* q^2 written with SP[], then reduced with the kinematic rules *)
sqExpand[q_] := Module[{qq = Expand[q], basis, cs},
  basis = Union[Cases[qq, _p | _l, {0, Infinity}]];
  cs = Coefficient[qq, #] & /@ basis;
  Expand[Sum[cs[[i]] cs[[j]] SP[basis[[i]], basis[[j]]], {i, Length[basis]}, {j, Length[basis]}]]];
ksq[q_, kr_] := With[{x = Expand[sqExpand[q] /. kr["Rules"]], sm = kr["Sum"]},
   If[sm === None, x,
    First[MinimalBy[{x, Expand[x /. u -> sm - s - t], Expand[x /. t -> sm - s - u], Expand[x /. s -> sm - t - u]},
      LeafCount]]]];


(* ::Section:: *)
(*Which vertices can appear*)


weakCompositions[v_, n_] := Join @@ (Permutations /@ (IntegerPartitions[v + n, {n}] - 1));

vertexSets[th_, exts_, loopsOpt_, orderOpt_] :=
 Module[{terms = th["Terms"], nT, degs, nE = Length[exts], loops, vlist, crules = None, vmax,
   nS = Length[th["SlotTypes"]], extCount, res = {}, expo, balancedQ},
  nT = Length[terms]; degs = terms[[All, "Degree"]];
  loops = Which[loopsOpt === Automatic && orderOpt === Automatic, {0},
    loopsOpt === Automatic, All, True, Flatten[{loopsOpt}]];
  (* coupling powers: couplings that are not listed must appear to power 0 *)
  If[MatchQ[orderOpt, {(_ -> _Integer) ..}],
   crules = Join[orderOpt, (# -> 0) & /@ Complement[
       Union[Flatten[Cases[#Coefficient, c_Symbol /; ! NumericQ[c], {0, Infinity}] & /@ terms]],
       orderOpt[[All, 1]]]]];
  If[IntegerQ[orderOpt] || VectorQ[orderOpt, IntegerQ],
   vlist = Flatten[{orderOpt}],
   vmax = Infinity;
   If[loops =!= All && Min[degs] >= 3, vmax = Min[vmax, 2 Max[loops] - 2 + nE]];
   If[crules =!= None && AllTrue[terms, Function[tm, Total[Exponent[tm["Coefficient"], #] & /@ crules[[All, 1]]] >= 1]],
    vmax = Min[vmax, Total[crules[[All, 2]]]]];
   If[vmax === Infinity, Message[FeynmanDiagrams::order]; Throw[$Failed, $fsTag]];
   vlist = Range[1, vmax]];
  vlist = Select[vlist, # >= 1 &];
  extCount = Table[Count[exts[[All, "Slot"]], sl], {sl, nS}];
  expo = If[crules === None, None, Table[Exponent[tm["Coefficient"], r[[1]]], {r, crules}, {tm, terms}]];
  balancedQ[free_] := AllTrue[Range[nS], Function[sl,
     If[th["SlotTypes"][[sl, 2]] == 0, EvenQ[free[[sl]]], free[[sl]] == free[[th["Conj"][[sl]]]]]]];
  Do[
   Do[
    Module[{free, nl},
     free = Total[c terms[[All, "Exponents"]]] - extCount;
     nl = (c . degs - nE)/2 - v + 1;
     If[AllTrue[free, NonNegative] && balancedQ[free] && IntegerQ[nl] && nl >= 0 &&
       (loops === All || MemberQ[loops, nl]) && (expo === None || expo . c == crules[[All, 2]]),
      AppendTo[res, Flatten[MapThread[ConstantArray, {Range[nT], c}]]]]],
    {c, weakCompositions[v, nT]}],
   {v, vlist}];
  res]


(* ::Section:: *)
(*Graph generation*)


(* Multigraphs with the given vertices. Edges are {tail, head, field}; for real fields tail <= head,
   for complex fields the arrow (particle flow) runs from the chi^* end (tail) to the chi end (head).
   Symmetry breaking: a free slot is paired with a vertex (not a particular leg), and among identical
   untouched vertices only the first one is tried. Duplicates are removed later by canonical labelling. *)

genGraphs[vt_List, slots_List, extSlots_List, conj_List, st_List] :=
 Module[{ctx},
  ctx = <|"V" -> Length[vt], "vt" -> vt, "ext" -> extSlots, "E" -> Length[extSlots],
    "conj" -> conj, "st" -> st, "nS" -> Length[st]|>;
  Flatten[Reap[attachExt[ctx, 1, slots, ConstantArray[False, Length[vt]], {}]][[2]], 1]]

breakSym[cands_List, touched_List, vt_List] :=
 Module[{seen = {}},
  Select[cands, Function[v,
    If[touched[[v]], True, If[MemberQ[seen, vt[[v]]], False, AppendTo[seen, vt[[v]]]; True]]]]]

attachExt[ctx_, j_, free_, touched_, att_] :=
 If[j > ctx["E"],
  pairInternal[ctx, free, touched, att, {}],
  With[{sl = ctx["ext"][[j]]},
   Do[
    attachExt[ctx, j + 1, ReplacePart[free, {v, sl} -> free[[v, sl]] - 1],
     ReplacePart[touched, v -> True], Append[att, v]],
    {v, breakSym[Select[Range[ctx["V"]], free[[#, sl]] > 0 &], touched, ctx["vt"]]}]]]

pairInternal[ctx_, free_, touched_, att_, edges_] :=
 Module[{u, tt, tc, cands},
  u = SelectFirst[Range[ctx["V"]], Total[free[[#]]] > 0 &, 0];
  If[u == 0,
   Sow[{att, Sort[edges]}],
   tt = SelectFirst[Range[ctx["nS"]], free[[u, #]] > 0 &];
   tc = ctx["conj"][[tt]];
   cands = Join[
     If[If[tc == tt, free[[u, tt]] >= 2, free[[u, tc]] >= 1], {u}, {}],
     breakSym[Select[Range[u + 1, ctx["V"]], free[[#, tc]] > 0 &], touched, ctx["vt"]]];
   Do[
    Module[{nf = free},
     nf[[u, tt]] -= 1; nf[[w, tc]] -= 1;
     pairInternal[ctx, nf, ReplacePart[touched, {u -> True, w -> True}], att,
      Append[edges, mkEdge[ctx, u, tt, w]]]],
    {w, cands}]]]

mkEdge[ctx_, u_, tt_, w_] :=
 With[{k = ctx["st"][[tt, 1]], sg = ctx["st"][[tt, 2]]},
  Switch[sg, 0, {Min[u, w], Max[u, w], k}, 1, {w, u, k}, -1, {u, w, k}]]

components[v_, pairs_] := ConnectedComponents[Graph[Range[v], UndirectedEdge @@@ pairs]];
connectedQ[v_, edges_] := If[edges === {}, v == 1, Length[components[v, edges[[All, 1 ;; 2]]]] == 1];

(* For every bridge (an internal line whose removal disconnects the diagram):
   the number of external legs on each side, sorted *)
bridges[v_, att_, edges_] :=
 If[edges === {}, {},
  Module[{pairs = edges[[All, 1 ;; 2]]},
   Flatten[Reap[Do[
       If[pairs[[i, 1]] != pairs[[i, 2]],
        With[{cc = components[v, Delete[pairs, i]]},
         If[Length[cc] == 2,
          Sow[Sort[{Count[att, _?(MemberQ[cc[[1]], #] &)], Count[att, _?(MemberQ[cc[[2]], #] &)]}]]]]],
       {i, Length[pairs]}]][[2]], 1]]]

keepQ[type_, tadpoles_, nE_, br_] :=
 And[
  Switch[type,
   "1PI", br === {},
   "Amputated", nE < 2 || FreeQ[br, {1, nE - 1}],   (* no self-energy on an external leg *)
   "Connected", True],
  tadpoles || FreeQ[br, {0, _}]]


(* ::Section:: *)
(*Canonical form and symmetry factor*)


rankList[list_List] := With[{a = AssociationThread[Union[list] -> Range[Length[Union[list]]]]}, a /@ list];

incSig[v_, col_, edges_, realQ_] :=
 Join @@ (Function[e, With[{a = e[[1]], b = e[[2]], k = e[[3]]},
      Which[
       a == v && b == v, {{0, k, 0}},
       a == v, {{If[realQ[[k]], 1, 2], k, col[[b]]}},
       b == v, {{If[realQ[[k]], 1, 3], k, col[[a]]}},
       True, {}]]] /@ edges)

codeFor[ord_, vt_, extAt_, edges_, realQ_] :=
 Module[{newIdx = ConstantArray[0, Length[ord]]},
  newIdx[[ord]] = Range[Length[ord]];
  {Transpose[{vt[[ord]], extAt[[ord]]}],
   Sort[Function[e, With[{a = newIdx[[e[[1]]]], b = newIdx[[e[[2]]]]},
        If[realQ[[e[[3]]]], {Min[a, b], Max[a, b], e[[3]]}, {a, b, e[[3]]}]]] /@ edges]}]

(* Colour refinement, then all orderings inside the colour classes. The smallest code is canonical;
   the number of orderings that reproduce it is the number of vertex automorphisms
   (external legs are labelled, so they stay fixed). *)
canon[vt_List, att_List, edges_List, realQ_List] :=
 Module[{nv = Length[vt], extAt, col, newcol, classes, ords, codes, best, pos},
  extAt = Table[Flatten[Position[att, v]], {v, nv}];
  col = rankList[Transpose[{vt, extAt}]];
  While[True,
   newcol = rankList[Table[{col[[v]], Sort[incSig[v, col, edges, realQ]]}, {v, nv}]];
   If[CountDistinct[newcol] == CountDistinct[col], Break[]];
   col = newcol];
  classes = GatherBy[SortBy[Range[nv], col[[#]] &], col[[#]] &];
  ords = Flatten /@ Tuples[Permutations /@ classes];
  codes = codeFor[#, vt, extAt, edges, realQ] & /@ ords;
  best = First[Sort[codes]];
  pos = Flatten[Position[codes, best, {1}, Heads -> False]];
  <|"Code" -> best, "AutV" -> Length[pos], "Order" -> ords[[First[pos]]]|>]

(* Automorphisms that fix every vertex: permutations of parallel lines (m!) and,
   for real fields, flips of self-loops (2 each) *)
multFactor[edges_, realQ_] :=
  Times @@ (If[#[[1, 1]] == #[[1, 2]] && realQ[[#[[1, 3]]]], 2^#[[2]] #[[2]]!, #[[2]]!] & /@ Tally[edges]);

graphsFor[th_, exts_, vt_] :=
 Module[{nv = Length[vt], raw},
  raw = genGraphs[vt, th["Terms"][[#, "Exponents"]] & /@ vt, exts[[All, "Slot"]], th["Conj"], th["SlotTypes"]];
  raw = Select[raw, connectedQ[nv, #[[2]]] &];
  DeleteDuplicatesBy[
   Function[g, Module[{c = canon[vt, g[[1]], g[[2]], th["RealQ"]], ord, newIdx},
      ord = c["Order"]; newIdx = ConstantArray[0, nv]; newIdx[[ord]] = Range[nv];
      <|"Code" -> c["Code"], "AutV" -> c["AutV"], "vt" -> vt[[ord]], "att" -> newIdx[[g[[1]]]],
        "edges" -> c["Code"][[2]]|>]] /@ raw,
   #Code &]]


(* ::Section:: *)
(*Momentum routing*)


(* Spanning tree from the vertex of the last external leg; each remaining line (chord) carries a
   loop momentum l[n] from tail to head; tree lines follow from momentum conservation.
   Returns the momentum of every line, flowing tail -> head. *)
route[nv_, att_, edges_, inflow_] :=
 Module[{nEd = Length[edges], root, adj, visited, parent, pedge, treeQ, order = {}, queue, chords, q, sub},
  root = att[[-1]];
  adj = ConstantArray[{}, nv];
  Do[With[{a = edges[[i, 1]], b = edges[[i, 2]]},
    If[a != b, AppendTo[adj[[a]], {i, b}]; AppendTo[adj[[b]], {i, a}]]], {i, nEd}];
  visited = ConstantArray[False, nv]; parent = ConstantArray[0, nv];
  pedge = ConstantArray[0, nv]; treeQ = ConstantArray[False, nEd];
  visited[[root]] = True; queue = {root};
  While[queue =!= {},
   With[{v = First[queue]},
    queue = Rest[queue]; AppendTo[order, v];
    Do[If[! visited[[iw[[2]]]],
      visited[[iw[[2]]]] = True; parent[[iw[[2]]]] = v; pedge[[iw[[2]]]] = iw[[1]];
      treeQ[[iw[[1]]]] = True; AppendTo[queue, iw[[2]]]], {iw, adj[[v]]}]]];
  chords = Select[Range[nEd], ! treeQ[[#]] &];
  q = ConstantArray[0, nEd];
  Do[q[[chords[[n]]]] = l[n], {n, Length[chords]}];
  (* momentum flowing into each vertex from external legs and chords *)
  sub = Table[Total[inflow[[Flatten[Position[att, v]]]]] +
     Sum[If[edges[[c, 2]] == v, q[[c]], 0] - If[edges[[c, 1]] == v, q[[c]], 0], {c, chords}], {v, nv}];
  Do[If[w != root,
    sub[[parent[[w]]]] += sub[[w]];
    q[[pedge[[w]]]] = If[edges[[pedge[[w]], 1]] == w, sub[[w]], -sub[[w]]]],
   {w, Reverse[order]}];
  Expand[q]]


(* ::Section:: *)
(*Amplitude of one diagram*)


chanName[x_] := Replace[x, {s -> "s", t -> "t", u -> "u", _ -> None}];

(* q^2 - m^2 with the sign of q chosen so that its first loop momentum enters with +1 *)
den[q_, m_] := With[{c = Coefficient[q, First[Union[Cases[q, _l, {0, Infinity}]]]]},
   Den[Expand[If[c < 0, -q, q]], m]];

(* one-loop integral over loop momentum lm of 1/Product[(q_i^2 - m_i^2)], in units of i/(16 pi^2) *)
oneLoopPV[lm_, qs_, ms_, eds_, kr_] :=
 Module[{n = Length[qs], c, r, from, to, ord, nxt, rr, mm, val, nm},
  c = Coefficient[qs, lm];
  r = Expand[c qs - lm];                 (* q_i = c_i (l + r_i),  c_i = +-1 *)
  from = MapThread[If[#1 > 0, #2[[1]], #2[[2]]] &, {c, eds}];
  to = MapThread[If[#1 > 0, #2[[2]], #2[[1]]] &, {c, eds}];
  ord = {1};                             (* walk around the loop in the direction of l *)
  Do[nxt = SelectFirst[Complement[Range[n], ord], from[[#]] == to[[Last[ord]]] &, Missing[]];
   If[MissingQ[nxt], ord = Range[n]; Break[]]; AppendTo[ord, nxt], {n - 1}];
  rr = r[[ord]]; mm = ms[[ord]]^2;
  val = Switch[n,
    1, A0[mm[[1]]],
    2, B0[ksq[rr[[2]] - rr[[1]], kr], Sequence @@ Sort[mm]],
    3, C0[ksq[rr[[2]] - rr[[1]], kr], ksq[rr[[3]] - rr[[2]], kr], ksq[rr[[3]] - rr[[1]], kr], Sequence @@ mm],
    4, D0[ksq[rr[[2]] - rr[[1]], kr], ksq[rr[[3]] - rr[[2]], kr], ksq[rr[[4]] - rr[[3]], kr],
     ksq[rr[[4]] - rr[[1]], kr], ksq[rr[[3]] - rr[[1]], kr], ksq[rr[[4]] - rr[[2]], kr], Sequence @@ mm],
    _, ScalarNPoint[rr, mm]];
  nm = Switch[n, 1, "tadpole", 2, "bubble" <> If[chanName[val[[1]]] =!= None, "(" <> chanName[val[[1]]] <> ")", ""],
    3, "triangle", 4, "box", 5, "pentagon", _, ToString[n] <> "-gon"];
  {val, nm}]

(* integrals over the loop momenta; factorises into independent loops where possible *)
loopIntegrals[{}, ___] := {1, {}};
loopIntegrals[qs_List, ms_List, eds_List, kr_] :=
 Module[{lofs, lms, grp, vals = {}, names = {}},
  lofs = Union[Cases[#, _l, {0, Infinity}]] & /@ qs;
  lms = Union[Flatten[lofs]];
  grp = ConnectedComponents[Graph[lms, Flatten[UndirectedEdge @@@ Partition[#, 2, 1] & /@ lofs]]];
  Do[
   With[{idx = Select[Range[Length[qs]], Intersection[lofs[[#]], g] =!= {} &]},
    If[Length[g] == 1,
     With[{r = oneLoopPV[g[[1]], qs[[idx]], ms[[idx]], eds[[idx]], kr]},
      AppendTo[vals, I/(16 Pi^2) r[[1]]]; AppendTo[names, r[[2]]]],
     AppendTo[vals, LoopIntegral[Sort[g], MapThread[den, {qs[[idx]], ms[[idx]]}]]];
     AppendTo[names, ToString[Length[g]] <> "-loop"]]],
   {g, SortBy[grp, Min[#[[All, 1]]] &]}];
  {Times @@ vals, names}]

amplitude[th_, vt_, edges_, q_, sym_, kr_, nE_] :=
 Module[{nEd = Length[edges], vf, ms, loopIdx, treeIdx, treeInv, treePart, lv, names, iM, chan, label, integrand},
  vf = Times @@ (th["Terms"][[#, "Factor"]] & /@ vt);
  ms = th["Fields"][[#[[3]], "Mass"]] & /@ edges;
  loopIdx = Select[Range[nEd], ! FreeQ[q[[#]], _l] &];
  treeIdx = Complement[Range[nEd], loopIdx];
  treeInv = ksq[q[[#]], kr] & /@ treeIdx;
  (* an on-shell propagator (self-energy on an external leg) is kept as I/Den[q, m] *)
  treePart = Times @@ MapThread[If[PossibleZeroQ[#1 - #2^2], I/Den[#3, #2], I/(#1 - #2^2)] &,
     {treeInv, ms[[treeIdx]], q[[treeIdx]]}];
  {lv, names} = loopIntegrals[q[[loopIdx]], ms[[loopIdx]], edges[[loopIdx]], kr];
  iM = vf/sym treePart I^Length[loopIdx] lv;
  integrand = vf/sym treePart Times @@ MapThread[I/den[#1, #2] &, {q[[loopIdx]], ms[[loopIdx]]}];
  chan = If[nE == 4, DeleteCases[DeleteDuplicates[chanName /@ treeInv], None], {}];
  label = Which[
    loopIdx === {} && treeIdx === {}, "contact",
    loopIdx === {}, If[Length[treeIdx] == 1 && chan =!= {}, chan[[1]] <> "-channel",
      "tree" <> If[chan =!= {}, " (" <> StringRiffle[chan, ","] <> ")", ""]],
    True, StringRiffle[names, " \[Times] "] <>
     If[chan =!= {}, " on " <> StringRiffle[chan, ","] <> "-channel", ""]];
  <|"VertexFactor" -> vf, "iM" -> Simplify[iM], "M" -> Simplify[-I iM], "Integrand" -> integrand,
    "Label" -> label, "LoopPropagators" -> Length[loopIdx], "TreeInvariants" -> treeInv|>]


(* ::Section:: *)
(*FeynmanDiagrams*)


Options[FeynmanDiagrams] = {"Loops" -> Automatic, "Order" -> Automatic, "Diagrams" -> "Amputated",
   "Tadpoles" -> True, "OnShell" -> Automatic};

FeynmanDiagrams[ScalarTheory[th_Association], proc_Rule, OptionsPattern[]] :=
 Catch[
  Module[{exts, nE, type, tad, onShell, kr, inflow, vsets, graphs, diags},
   exts = parseProcess[th, proc];
   nE = Length[exts];
   type = OptionValue["Diagrams"]; tad = TrueQ[OptionValue["Tadpoles"]];
   If[! MemberQ[{"Amputated", "1PI", "Connected"}, type], Message[FeynmanDiagrams::dtype, type]; Throw[$Failed, $fsTag]];
   onShell = Replace[OptionValue["OnShell"], Automatic -> (nE =!= 2 && type =!= "Connected")];
   kr = kinRules[exts, onShell];
   (* momentum flowing into the diagram on each leg; a 2-point function has one momentum p1 *)
   inflow = Which[nE == 1, {0}, nE == 2, {p[1], -p[1]},
     True, (If[#Direction == "In", 1, -1] #Momentum) & /@ exts];
   exts = MapThread[Append[#1, "Inflow" -> #2] &, {exts, inflow}];
   If[nE == 2, exts = Append[#, "Momentum" -> p[1]] & /@ exts];
   vsets = vertexSets[th, exts, OptionValue["Loops"], OptionValue["Order"]];
   graphs = DeleteDuplicatesBy[Join @@ (graphsFor[th, exts, #] & /@ vsets), #Code &];
   graphs = Select[graphs, keepQ[type, tad, nE, bridges[Length[#vt], #att, #edges]] &];
   diags = Function[g, Module[{q, sym, amp},
       q = route[Length[g["vt"]], g["att"], g["edges"], exts[[All, "Inflow"]]];
       sym = g["AutV"] multFactor[g["edges"], th["RealQ"]];
       amp = amplitude[th, g["vt"], g["edges"], q, sym, kr, nE];
       Join[<|"Theory" -> th, "Process" -> proc, "ProcessTeX" -> processTeX[exts],
         "Externals" -> MapThread[Append[#1, "Vertex" -> #2] &, {exts, g["att"]}],
         "VertexTypes" -> g["vt"], "Edges" -> g["edges"], "Momenta" -> q,
         "Loops" -> Length[Union[Cases[q, _l, {0, Infinity}]]], "SymmetryFactor" -> sym,
         "Code" -> g["Code"]|>, amp]]] /@ graphs;
   diags = SortBy[diags, {#Loops, Length[#VertexTypes], Sort[#VertexTypes], -#LoopPropagators, #Code} &];
   diags = MapIndexed[FeynmanDiagram[Prepend[#1, "Id" -> First[#2]]] &, diags];
   If[diags === {}, Message[FeynmanDiagrams::none]];
   diags],
  $fsTag]


(* ::Section:: *)
(*Tables and totals*)


DiagramTable[d_FeynmanDiagram, rest___] := DiagramTable[{d}, rest];
DiagramTable[diags : {__FeynmanDiagram}, key_String : "iM"] :=
 Grid[Prepend[
   Function[a, {a["Id"], a["Label"], a["Loops"], a["SymmetryFactor"], TraditionalForm[a[key]]}] /@ diags[[All, 1]],
   Style[#, Bold] & /@ {"#", "diagram", "loops", "S", key}],
  Frame -> All, Alignment -> {Left, Center}, Spacings -> {1.5, 0.8},
  Background -> {None, {GrayLevel[0.92]}}]

Amplitude[diags : {___FeynmanDiagram}, form_String : "M"] :=
 Collect[Total[#[[1, If[form === "iM", "iM", "M"]]] & /@ diags],
  {_A0, _B0, _C0, _D0, _ScalarNPoint, _LoopIntegral}, Simplify]


(* ::Section:: *)
(*Loop functions*)


PVExplicit[expr_] :=
 Module[{x = \[FormalX], y = \[FormalY], z = \[FormalZ], eps = \[FormalEpsilon]},
  expr /. {
    A0[0] -> 0,
    A0[m2_] :> m2 (DeltaUV + 1 - Log[m2/MuR^2]),
    B0[p2_, m12_, m22_] :> DeltaUV -
      Inactive[Integrate][Log[(x m12 + (1 - x) m22 - x (1 - x) p2 - I eps)/MuR^2], {x, 0, 1}],
    C0[p1_, p2_, p3_, m1_, m2_, m3_] :> -Inactive[Integrate][
       1/(x m1 + y m2 + (1 - x - y) m3 - x y p1 - y (1 - x - y) p2 - x (1 - x - y) p3 - I eps),
       {x, 0, 1}, {y, 0, 1 - x}],
    D0[p1_, p2_, p3_, p4_, s12_, s23_, m1_, m2_, m3_, m4_] :> With[{w = 1 - x - y - z},
      Inactive[Integrate][
       1/(x m1 + y m2 + z m3 + w m4 - x y p1 - y z p2 - z w p3 - x w p4 - x z s12 - y w s23 - I eps)^2,
       {x, 0, 1}, {y, 0, 1 - x}, {z, 0, 1 - x - y}]]}]

Options[EvaluateLoops] = {"DeltaUV" -> Automatic, "Mu" -> Automatic};
EvaluateLoops[expr_, rules : {___Rule} : {}, OptionsPattern[]] :=
 Module[{e, fix},
  fix[x_] := Module[{y = x},
    If[OptionValue["Mu"] =!= Automatic, y = y /. MuR -> OptionValue["Mu"]];
    If[OptionValue["DeltaUV"] =!= Automatic, y = y /. DeltaUV -> OptionValue["DeltaUV"]]; y];
  e = fix[expr /. rules];
  e = e /. {
     A0[a_?NumericQ] :> nA0[a],
     B0[a_?NumericQ, b_?NumericQ, c_?NumericQ] :> nB0[a, b, c],
     C0[x__?NumericQ] /; Length[{x}] == 6 :> nC0[N[{x}]],
     D0[x__?NumericQ] /; Length[{x}] == 10 :> nD0[N[{x}]]};
  Chop[N[fix[e]]]]

nA0[m2_] := If[m2 == 0, 0, m2 (DeltaUV + 1 - Log[N[m2]] + Log[MuR^2])];

(* B0 = DeltaUV - Integral_0^1 dx Log[(D(x) - i0)/mu^2],  D(x) = x m1 + (1-x) m2 - x(1-x) p2 *)
nB0[p2_, m1_, m2_] := If[p2 == 0 && m1 == 0 && m2 == 0, 0, DeltaUV + Log[MuR^2] - bInt[N[p2], N[m1], N[m2]]];
bInt[p2_, m1_, m2_] := Module[{x, dd, sol, roots},
  dd = x m1 + (1 - x) m2 - x (1 - x) p2;
  sol = Quiet[NSolve[dd == 0, x, Reals]];
  roots = If[ListQ[sol] && sol =!= {}, Select[x /. sol, NumericQ[#] && 0 < # < 1 &], {}];
  NIntegrate[Log[Abs[dd]] - I Pi UnitStep[-dd], Evaluate[Join[{x, 0}, Sort[roots], {1}]],
   Method -> "GlobalAdaptive", MaxRecursion -> 30]]

(* Feynman-parameter integral of jac/F^k over the unit cube. If F changes sign (above a threshold)
   the contour is deformed, x_j -> x_j - i lam x_j (1 - x_j) dF/dx_j, which gives Im F < 0 as the
   -i eps prescription requires. *)
cubeIntegral[F_, jac_, k_, vars_, lam0_] :=
 Module[{grid, vals, grad, gmax, lam, zs, G, J},
  grid = Tuples[Range[0, 10]/10., Length[vars]];
  vals = (F /. Thread[vars -> #]) & /@ grid;
  If[Min[vals] > 0 || Max[vals] < 0,
   NIntegrate[Evaluate[jac/F^k], Evaluate[Sequence @@ ({#, 0, 1} & /@ vars)], MaxRecursion -> 20],
   grad = D[F, {vars}];
   gmax = Max[Abs[Flatten[(grad /. Thread[vars -> #]) & /@ grid]]];
   lam = lam0/Max[gmax, 10^-12];
   zs = vars - I lam vars (1 - vars) grad;
   G = (jac/F^k) /. Thread[vars -> zs];
   J = Det[D[zs, {vars}]];
   NIntegrate[Evaluate[G J], Evaluate[Sequence @@ ({#, 0, 1} & /@ vars)], MaxRecursion -> 20,
    PrecisionGoal -> 6, Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 10^4}]]]

(* simplex -> unit cube: x = X, y = (1-X) Y, z = (1-X)(1-Y) Z *)
nC0[{p1_, p2_, p3_, m1_, m2_, m3_}] :=
 Module[{X, Y, x, y},
  x = X; y = (1 - X) Y;
  -cubeIntegral[Expand[x m1 + y m2 + (1 - x - y) m3 - x y p1 - y (1 - x - y) p2 - x (1 - x - y) p3],
    1 - X, 1, {X, Y}, 0.5]]

nD0[{p1_, p2_, p3_, p4_, s12_, s23_, m1_, m2_, m3_, m4_}] :=
 Module[{X, Y, Z, x, y, z, w},
  x = X; y = (1 - X) Y; z = (1 - X) (1 - Y) Z; w = 1 - x - y - z;
  cubeIntegral[Expand[x m1 + y m2 + z m3 + w m4 - x y p1 - y z p2 - z w p3 - x w p4 - x z s12 - y w s23],
   (1 - X)^2 (1 - Y), 2, {X, Y, Z}, 1.]]


(* ::Section:: *)
(*Export and drawing*)


texStr[x_] := ToString[TeXForm[x]];
(* lines of real fields carry no arrow, so their momentum label is shown with a leading + sign *)
posMom[q_] := With[{c = Cases[Expand[q], _p | _l, {0, Infinity}]},
   If[c =!= {} && Coefficient[Expand[q], First[SortBy[c, {Head[#] =!= l &, First}]]] < 0, -q, q]];
texAmp[x_] := texStr[x /. LoopIntegral[ls_, _] :> Subscript[\[ScriptCapitalI], Length[ls]]];

fieldJSON[f_] := <|"index" -> f["Index"], "name" -> ToString[f["Symbol"]], "tex" -> texStr[f["Symbol"]],
   "type" -> ToLowerCase[f["Type"]], "mass_tex" -> texStr[f["Mass"]]|>;
termJSON[tm_] := <|"index" -> tm["Index"], "factor_tex" -> texStr[tm["Factor"]],
   "term_tex" -> texStr[tm["Coefficient"] tm["Monomial"]]|>;
extJSON[e_] := <|"id" -> e["Index"], "field" -> e["Field"], "anti" -> e["Anti"],
   "direction" -> ToLowerCase[e["Direction"]], "vertex" -> e["Vertex"], "momentum_tex" -> texStr[e["Momentum"]],
   "arrow" -> Switch[e["SlotType"][[2]], 0, "none", 1, "to_vertex", -1, "from_vertex"]|>;
diagJSON[a_] := <|"id" -> a["Id"], "label" -> a["Label"], "symmetry_factor" -> a["SymmetryFactor"],
   "loops" -> a["Loops"], "process_tex" -> a["ProcessTeX"],
   "vertices" -> MapIndexed[<|"id" -> First[#2], "term" -> #1|> &, a["VertexTypes"]],
   "externals" -> (extJSON /@ a["Externals"]),
   "edges" -> MapThread[<|"tail" -> #1[[1]], "head" -> #1[[2]], "field" -> #1[[3]],
       "directed" -> (a["Theory"]["Fields"][[#1[[3]], "Type"]] == "Complex"),
       "momentum_tex" -> texStr[If[a["Theory"]["Fields"][[#1[[3]], "Type"]] == "Real", posMom[#2], #2]]|> &,
     {a["Edges"], a["Momenta"]}],
   "amplitude_tex" -> texAmp[a["iM"]]|>;

ExportDiagrams[file_String, d_FeynmanDiagram] := ExportDiagrams[file, {d}];
ExportDiagrams[file_String, diags : {__FeynmanDiagram}] :=
 Module[{a1 = diags[[1, 1]], th},
  th = a1["Theory"];
  Export[file, <|"title" -> a1["ProcessTeX"], "fields" -> (fieldJSON /@ th["Fields"]),
    "terms" -> (termJSON /@ th["Terms"]), "diagrams" -> (diagJSON /@ diags[[All, 1]])|>, "RawJSON"]]

findPython[] := findPython[] = Module[{cands},
   cands = Join[
     Select[{"/opt/anaconda3/bin/python3", FileNameJoin[{$HomeDirectory, "anaconda3", "bin", "python3"}],
       FileNameJoin[{$HomeDirectory, "miniconda3", "bin", "python3"}], "/opt/miniconda3/bin/python3",
       "/opt/homebrew/bin/python3", "/usr/local/bin/python3", "/usr/bin/python3"}, FileExistsQ],
     {"python3", "python"}];
   SelectFirst[cands, Quiet[Check[RunProcess[{#, "-c", "import numpy, matplotlib"}, "ExitCode"], 1]] === 0 &, $Failed]];

Options[DrawDiagrams] = {"File" -> None, "Columns" -> Automatic, "ShowMomenta" -> False,
   "ShowAmplitude" -> False, "Separate" -> None, "LineStyle" -> Automatic};
DrawDiagrams[d_FeynmanDiagram, opts : OptionsPattern[]] := DrawDiagrams[{d}, opts];
DrawDiagrams[diags : {__FeynmanDiagram}, OptionsPattern[]] :=
 Module[{py, json, png, args, res, file = OptionValue["File"], sep = OptionValue["Separate"]},
  py = If[StringQ[$FeynScalarPython], $FeynScalarPython, findPython[]];
  If[! StringQ[py], Message[DrawDiagrams::nopy]; Return[$Failed]];
  json = FileNameJoin[{$TemporaryDirectory, "feynscalar-" <> CreateUUID[] <> ".json"}];
  png = StringDrop[json, -5] <> ".png";
  ExportDiagrams[json, diags];
  args = {py, FileNameJoin[{$FeynScalarDirectory, "draw_feynman.py"}], json, "-o", png};
  If[StringQ[file], args = Join[args, {"-o", ExpandFileName[file]}]];
  If[StringQ[sep], args = Join[args, {"--separate", ExpandFileName[sep]}]];
  If[IntegerQ[OptionValue["Columns"]], args = Join[args, {"--cols", ToString[OptionValue["Columns"]]}]];
  If[TrueQ[OptionValue["ShowMomenta"]], AppendTo[args, "--momenta"]];
  If[TrueQ[OptionValue["ShowAmplitude"]], AppendTo[args, "--amplitude"]];
  If[StringQ[OptionValue["LineStyle"]], args = Join[args, {"--real-style", ToLowerCase[OptionValue["LineStyle"]]}]];
  res = RunProcess[args];
  If[res["ExitCode"] =!= 0, Message[DrawDiagrams::pyerr, res["StandardError"]]; Return[$Failed]];
  Import[png]]

End[];
EndPackage[];
