#!/usr/bin/env wolframscript
(* FeynScalar examples and checks against textbook results.
   Run:  wolframscript -file examples.wls
   Pictures are written to ./examples/                                                    *)

dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "FeynScalar.wl"}]];
fig[name_] := FileNameJoin[{dir, "examples", name}];
Quiet[CreateDirectory[FileNameJoin[{dir, "examples"}]]];

pass = 0; fail = 0;
check[name_, value_, expected_] :=
  If[value === expected || (! ListQ[value] && ! AssociationQ[value] && Simplify[value - expected] === 0),
   pass++; Print["  ok    ", name],
   fail++; Print["  FAIL  ", name, "\n        got      ", InputForm[value], "\n        expected ", InputForm[expected]]];
show[label_, diags_] := (Print["\n", label, "  (", Length[diags], " diagrams)"];
   Scan[Print["   #", #["Id"], "  ", #["Label"], "   S = ", #["SymmetryFactor"], "   iM = ", InputForm[#["iM"]]] &, diags]);

(* ------------------------------------------------------------------ 1. phi^4 *)
Print["\n=== real scalar, Lint = -lambda/4! phi^4 ==="];
phi4 = ScalarTheory[{RealScalar[\[Phi], m]}, -\[Lambda]/4! \[Phi]^4];

tree = FeynmanDiagrams[phi4, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}];
oneLoop = FeynmanDiagrams[phi4, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}, "Loops" -> 1];
show["phi phi -> phi phi, one loop", oneLoop];
check["tree M = -lambda", Amplitude[tree], -\[Lambda]];
check["one-loop M (Peskin 10.20)", Amplitude[oneLoop],
  \[Lambda]^2/(32 Pi^2) (B0[s, m^2, m^2] + B0[t, m^2, m^2] + B0[u, m^2, m^2])];

se = FeynmanDiagrams[phi4, {\[Phi]} -> {\[Phi]}, "Loops" -> {1, 2}];
show["self-energy to two loops", se];
check["self-energy symmetry factors 2, 4, 6", #["SymmetryFactor"] & /@ se, {2, 4, 6}];
check["one-loop self-energy -i Sigma", se[[1]]["iM"], I \[Lambda] A0[m^2]/(32 Pi^2)];

twoLoop = FeynmanDiagrams[phi4, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}, "Loops" -> 2];
check["two-loop 4-point: 12 diagrams, S = {2 x9, 4 x3}", Sort[#["SymmetryFactor"] & /@ twoLoop], Sort[Join[ConstantArray[2, 9], ConstantArray[4, 3]]]];

DrawDiagrams[Join[tree, oneLoop], "File" -> fig["phi4_2to2_1loop.png"], "ShowAmplitude" -> True];
DrawDiagrams[se, "File" -> fig["phi4_selfenergy.png"], "ShowMomenta" -> True];
DrawDiagrams[twoLoop, "File" -> fig["phi4_2to2_2loop.png"]];

(* numbers: MS-bar (DeltaUV = 0), mu = m = 1, lambda = 1, s = 5, t = -1, u = 4 - s - t *)
Print["  M(s=5, t=-1) one loop, MS-bar: ",
  EvaluateLoops[Amplitude[Join[tree, oneLoop]], {\[Lambda] -> 1, m -> 1, s -> 5, t -> -1, u -> 0}, "DeltaUV" -> 0, "Mu" -> 1]];

(* ------------------------------------------------------------------ 2. phi^3 *)
Print["\n=== real scalar, Lint = -g/3! phi^3 ==="];
phi3 = ScalarTheory[{RealScalar[\[Phi], m]}, -g/3! \[Phi]^3];
t3 = FeynmanDiagrams[phi3, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}];
check["tree M = -g^2 (1/(s-m^2) + 1/(t-m^2) + 1/(u-m^2))", Amplitude[t3],
  -g^2 (1/(s - m^2) + 1/(t - m^2) + 1/(u - m^2))];
l3 = FeynmanDiagrams[phi3, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}, "Loops" -> 1];
check["one loop: 3 boxes, 6 triangles, 3 bubbles, 3 tadpoles",
  Counts[First[StringSplit[#["Label"]]] & /@ l3], <|"box" -> 3, "triangle" -> 6, "bubble(s)" -> 1, "bubble(t)" -> 1, "bubble(u)" -> 1, "tadpole" -> 3|>];
check["without tadpoles: 12", Length[FeynmanDiagrams[phi3, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}, "Loops" -> 1, "Tadpoles" -> False]], 12];
check["1PI: 3 boxes", Length[FeynmanDiagrams[phi3, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}, "Loops" -> 1, "Diagrams" -> "1PI"]], 3];
DrawDiagrams[Join[t3, l3], "File" -> fig["phi3_2to2_1loop.png"]];

(* ------------------------------------------------------------------ 3. complex scalar *)
Print["\n=== complex scalar, Lint = -lambda |chi|^4 ==="];
cs = ScalarTheory[{ComplexScalar[\[Chi], M]}, -\[Lambda] Abs[\[Chi]]^4];
cc = FeynmanDiagrams[cs, {\[Chi], \[Chi]} -> {\[Chi], \[Chi]}, "Loops" -> {0, 1}];
show["chi chi -> chi chi", cc];
check["chi chi -> chi chi: S = 1, 2, 1, 1", #["SymmetryFactor"] & /@ cc, {1, 2, 1, 1}];
check["contact vertex -4 i lambda", cc[[1]]["iM"], -4 I \[Lambda]];
ca = FeynmanDiagrams[cs, {\[Chi], Anti[\[Chi]]} -> {\[Chi], Anti[\[Chi]]}, "Loops" -> {0, 1}];
check["crossing chi chibar: S(s) = 1, S(u) = 2", #["SymmetryFactor"] & /@ ca, {1, 1, 1, 2}];
DrawDiagrams[Join[cc, ca], "File" -> fig["complex_phi4.png"], "ShowMomenta" -> True];

(* ------------------------------------------------------------------ 4. scalar Yukawa *)
Print["\n=== scalar Yukawa, Lint = -g phi chi^* chi ==="];
yuk = ScalarTheory[{RealScalar[\[Phi], m], ComplexScalar[\[Chi], M, \[Chi]b]}, -g \[Phi] \[Chi]b \[Chi]];
ann = FeynmanDiagrams[yuk, {\[Chi], \[Chi]b} -> {\[Chi], \[Chi]b}];
check["chi chibar -> chi chibar: M = -g^2 (1/(s-m^2) + 1/(t-m^2))", Amplitude[ann], -g^2 (1/(s - m^2) + 1/(t - m^2))];
sca = FeynmanDiagrams[yuk, {\[Chi], \[Chi]} -> {\[Chi], \[Chi]}];
check["chi chi -> chi chi: M = -g^2 (1/(t-m^2) + 1/(u-m^2))", Amplitude[sca], -g^2 (1/(t - m^2) + 1/(u - m^2))];
pair = FeynmanDiagrams[yuk, {\[Chi], \[Chi]b} -> {\[Phi], \[Phi]}];
check["chi chibar -> phi phi: M = -g^2 (1/(t-M^2) + 1/(u-M^2))", Amplitude[pair], -g^2 (1/(t - M^2) + 1/(u - M^2))];
dec = FeynmanDiagrams[yuk, {\[Phi]} -> {\[Chi], \[Chi]b}];
check["phi -> chi chibar: M = -g", Amplitude[dec], -g];
seY = Join[FeynmanDiagrams[yuk, {\[Phi]} -> {\[Phi]}, "Loops" -> 1], FeynmanDiagrams[yuk, {\[Chi]} -> {\[Chi]}, "Loops" -> 1]];
show["self-energies", seY];
check["phi self-energy: charged loop, S = 1", seY[[1]]["iM"], I g^2 B0[SP[p[1], p[1]], M^2, M^2]/(16 Pi^2)];
DrawDiagrams[Join[ann, sca, pair, dec], "File" -> fig["yukawa_tree.png"]];
DrawDiagrams[seY, "File" -> fig["yukawa_selfenergy.png"], "ShowMomenta" -> True];

(* ------------------------------------------------------------------ 5. several couplings *)
Print["\n=== phi^3 + phi^4, order g^2 lambda ==="];
mix = ScalarTheory[{RealScalar[\[Phi], m]}, -g/3! \[Phi]^3 - \[Lambda]/4! \[Phi]^4];
gl = FeynmanDiagrams[mix, {\[Phi], \[Phi]} -> {\[Phi], \[Phi]}, "Order" -> {g -> 2, \[Lambda] -> 1}];
check["g^2 lambda: 33 diagrams, all one loop", {Length[gl], Union[#["Loops"] & /@ gl]}, {33, {1}}];

Print["\n", pass, " checks passed, ", fail, " failed.  Pictures in ", FileNameJoin[{dir, "examples"}]];
