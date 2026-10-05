(* ::Package:: *)

(* Tests for QNMRadialKN and Kerr-Newman QNMRadialFunction objects.
   Run from a fresh kernel with
     TestReport["/path/to/QNM-Mathematica/Tests/QNMRadialKN.wlt"]
   Runtime: about half a minute.

   The radial function is Hintz's w (arXiv:2609.33661, eqs. (14)-(15)), normalised to w(r+) = 1. It is evaluated from
   four pieces: the horizon series (r+ <= r <= r0), the horizon solution integrated to rmatch, the outgoing solution
   continued outwards to rfar, and the asymptotic series beyond rfar. The checks below use finite differences of the
   returned values, so that the residual of the radial system tests the values themselves (the package computes
   derivatives from the system). *)

$QNMRepository = ParentDirectory[DirectoryName[$TestFileName]];
PacletDirectoryLoad[$QNMRepository];
Get["QNM`"];

(* Hintz's radial system written out (M = 1):
     d/dr w = (R1 + I sg lam) ws,  Delta d/dr ws + ((1 - 3 sg)/2 Delta' + 2 I Kt) ws = (R1 - I sg lam) w,
   R1 = -sg 2 om Q (r - Q^2)/(1 - Q^2) + 3 I/(2 Q), Kt = -(r^2 + a^2) om + a m - xi0 Delta, xi0 = 2 om Q^2/(1 - Q^2).
   Relative residuals with 4th-order central differences of step h. *)
fd[f_, r_, h_] := (f[r - 2 h] - 8 f[r - h] + 8 f[r + h] - f[r + 2 h])/(12 h);
hintzResidual[R_, r_, h_ : 10^-3] := Module[{Q = R["Q"], a = R["a"], m = R["m"], om = R["\[Omega]"], lam = R["Eigenvalue"],
    sg = R["SpinSystem"], Rs = R["w\[Sharp]"], De, xi0, R1, Kt, c, w, ws, dw, dws},
  De = r^2 - 2 r + a^2 + Q^2;
  xi0 = 2 om Q^2/(1 - Q^2);
  R1 = -sg 2 om Q (r - Q^2)/(1 - Q^2) + 3 I/(2 Q);
  Kt = -(r^2 + a^2) om + a m - xi0 De;
  c = (1 - 3 sg)/2 (2 r - 2) + 2 I Kt;
  {w, ws} = {R[r], Rs[r]};
  {dw, dws} = {fd[R, r, h], fd[Rs, r, h]};
  Max[Abs[dw - (R1 + I sg lam) ws]/(Abs[dw] + Abs[(R1 + I sg lam) ws]),
    Abs[De dws + c ws - (R1 - I sg lam) w]/(Abs[De dws] + Abs[c ws] + Abs[(R1 - I sg lam) w])]];

(* piece boundaries r0, rmatch, rfar (internal data of the object) *)
pieces[R_] := Lookup[R[[1]]["RadialFunction"][[1]], {"r0", "rmatch", "rfar"}];
(* sample radii: two in each piece (the horizon series, the horizon solution, the outgoing solution, the asymptotic series) *)
samples[R_] := Module[{rp = R["Domain"][[1]], r0, rm, rf},
  {r0, rm, rf} = pieces[R];
  {rp + 0.4 (r0 - rp), rp + 0.9 (r0 - rp), r0 + 0.3 (rm - r0), rm - 0.5, rm + 2, 0.7 rf, 1.3 rf, 5 rf}];
(* relative jump across r = b, with the first-order change removed *)
jump[f_, b_, eps_ : 10^-7] := Abs[f[b + eps] - f[b - eps] - 2 eps f'[b]]/Abs[f[b]];

Rkn = QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5];
Rrn = QNMRadialKN[-1, 1, 0, 0, 0., 0.5];
Rnx = QNMRadialKN[-2, 2, 2, 0, N[0.70709513674232926], N[0.70709513674232926]];   (* 1 - a^2 - Q^2 = 3.3*10^-5 *)
all = {Rkn, Rrn, Rnx};

VerificationTest[
  StringStartsQ[FindFile["QNM`"], $QNMRepository]
  ,
  True
  ,
  TestID -> "KNRadial-LoadsRepositoryVersion"
]


(* ::Section:: *)
(*Object and keys*)


VerificationTest[
  {Head[#], SubsetQ[Keys[#], {"s", "l", "m", "n", "a", "Q", "\[Omega]", "Eigenvalue", "Family", "SpinSystem", "Component",
     "Method", "Amplitudes", "Coordinates", "Domain"}]} & /@ all
  ,
  ConstantArray[{QNMRadialFunction, True}, 3]
  ,
  TestID -> "KNRadial-keys"
]

VerificationTest[
  {#["\[Omega]"] - QNMFrequencyKN[#["s"], #["l"], #["m"], #["n"], #["a"], #["Q"]]} & /@ all
  ,
  {{0.}, {0.}, {0.}}
  ,
  SameTest -> (Max[Abs[#1 - #2]] < 10^-13 &),
  TestID -> "KNRadial-frequency-consistent-with-QNMFrequencyKN"
]

VerificationTest[
  {Rkn["Q"], Rkn["Family"], Rrn["Family"], Rkn["SpinSystem"], Rkn["Component"], Rkn["w\[Sharp]"]["Component"],
    Rkn["Coordinates"], Rkn["Method"], Rkn["RadialFunction"], Rkn["Eigenvalue"]}
  ,
  {0.5, "gravitational-led", "electromagnetic-led", 1, "w", "w\[Sharp]", "BL", "HintzSeparated",
    Missing["KeyAbsent", "RadialFunction"], 3.159184554208 + 0.067313088142 I}
  ,
  SameTest -> (Most[#1] === Most[#2] && Abs[Last[#1] - Last[#2]] < 10^-9 &),
  TestID -> "KNRadial-accessors"
]


(* ::Section:: *)
(*Normalisation, smoothness at the horizon, radial system*)


VerificationTest[
  {#[#["Domain"][[1]]], #["Amplitudes"]["\[ScriptCapitalH]"]} & /@ all
  ,
  {{1., 1.}, {1., 1.}, {1., 1.}}
  ,
  SameTest -> (Max[Abs[#1 - #2]] < 10^-13 &),
  TestID -> "KNRadial-normalisation-at-horizon"
]

(* w is smooth at r+: its derivatives at r+ are finite and continuous *)
VerificationTest[
  Max[Table[With[{rp = R["Domain"][[1]]}, Abs[R'[rp + 10^-8] - R'[rp]]/Abs[R'[rp]] + Abs[R''[rp + 10^-8] - R''[rp]]/Abs[R''[rp]]], {R, all}]] < 10^-5
  ,
  True
  ,
  TestID -> "KNRadial-smooth-at-horizon"
]

VerificationTest[
  Max[Table[hintzResidual[R, r], {R, {Rkn, Rrn}}, {r, samples[R]}]] < 10^-8
  ,
  True
  ,
  TestID -> "KNRadial-radial-system-residual"
]

VerificationTest[
  Max[Table[hintzResidual[Rnx, r, 10^-4], {r, Rest[samples[Rnx]]}]] < 10^-6
  ,
  True
  ,
  TestID -> "KNRadial-radial-system-residual-near-extremal"
]

VerificationTest[
  Max[Table[Max[jump[R, b], jump[R["w\[Sharp]"], b]], {R, all}, {b, pieces[R]}]] < 10^-8
  ,
  True
  ,
  TestID -> "KNRadial-continuity-across-pieces"
]


(* ::Section:: *)
(*Outgoing behaviour and coordinates*)


(* hyperboloidal function H = Exp[-I k (r - r+)] (r/r+)^-p w tends to the amplitude I like 1 + O(1/r) *)
VerificationTest[
  Table[Module[{H = QNMRadialKN[R["s"], R["l"], R["m"], R["n"], R["a"], R["Q"], "Frequency" -> R["\[Omega]"],
       Method -> {"HintzSeparated", "Coordinates" -> "Hyperboloidal"}], Ic = R["Amplitudes"]["\[ScriptCapitalI]"], d3, d6},
     d3 = Abs[H[10.^3]/Ic - 1]; d6 = Abs[H[10.^6]/Ic - 1];
     500 < d3/d6 < 2000 && d6 < 10^-4], {R, all}]
  ,
  {True, True, True}
  ,
  TestID -> "KNRadial-outgoing-at-infinity"
]

VerificationTest[
  Module[{k = Rkn["AsymptoticExponents"]["k"], p = Rkn["AsymptoticExponents"]["p"], rp = Rkn["Domain"][[1]],
     H = QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5, Method -> {"HintzSeparated", "Coordinates" -> "Hyperboloidal"}],
     S = QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5, Method -> {"HintzSeparated", "Coordinates" -> "CompactifiedHyperboloidal"}]},
    {p - 4 I Rkn["\[Omega]"], H[50.] - Exp[-I k (50 - rp)] (50/rp)^-p Rkn[50.], S[1/50.] - H[50.], S[0] - Rkn["Amplitudes"]["\[ScriptCapitalI]"],
     S'[0.02] - (-50.^2 H'[50.])}]
  ,
  {0, 0, 0, 0, 0}
  ,
  SameTest -> (Max[Abs[#1 - #2]] < 10^-9 &),
  TestID -> "KNRadial-coordinates"
]


(* ::Section:: *)
(*Derivatives and symmetries*)


VerificationTest[
  Max[Table[Abs[fd[Rkn, r, 10^-3] - Rkn'[r]]/Abs[Rkn'[r]] + Abs[fd[Rkn', r, 10^-3] - Rkn''[r]]/Abs[Rkn''[r]], {r, samples[Rkn]}]] < 10^-8
  ,
  True
  ,
  TestID -> "KNRadial-derivatives"
]

VerificationTest[
  Module[{Rm = QNMRadialKN[-2, 2, 2, 0, 0.6, -0.5]},
    {Rm[3.] - Rkn[3.], Rm["w\[Sharp]"][3.] + Rkn["w\[Sharp]"][3.], Rm["Eigenvalue"] + Rkn["Eigenvalue"], Rm["Q"] + 0.5}]
  ,
  {0, 0, 0, 0}
  ,
  SameTest -> (Max[Abs[#1 - #2]] < 10^-10 &),
  TestID -> "KNRadial-charge-reflection"
]

VerificationTest[
  Module[{R = QNMRadialKN[-1, 1, 0, 0, 0., 0.5, "Frequency" -> Rrn["\[Omega]"]]}, R[{2.5, 10., 100.}] - Rrn[{2.5, 10., 100.}]]
  ,
  {0, 0, 0}
  ,
  SameTest -> (Max[Abs[#1 - #2]] < 10^-12 &),
  TestID -> "KNRadial-frequency-option"
]


(* ::Section:: *)
(*Display and Kerr objects*)


VerificationTest[
  Module[{b = ToBoxes[Rkn]}, {Head[b], !FreeQ[b, "\"Q: \""], !FreeQ[b, "\"Family: \""]}]
  ,
  {InterpretationBox, True, True}
  ,
  TestID -> "KNRadial-summary-box"
]

VerificationTest[
  Module[{R = QNMRadial[-2, 2, 2, 0, 0.3]}, {FreeQ[ToBoxes[R], "\"Q: \""], R["Q"], R["w"], Head[R[3.]]}]
  ,
  {True, Missing["KeyAbsent", "Q"], Missing["KeyAbsent", "w"], Complex}
  ,
  TestID -> "KNRadial-Kerr-objects-unchanged"
]

VerificationTest[
  Rkn[1.]
  ,
  Indeterminate
  ,
  {QNMRadialFunction::dmval},
  TestID -> "KNRadial-domain"
]


(* ::Section:: *)
(*Small charge*)


(* 24-digit reference values computed with this package. For the electromagnetic-led family in the default spin
   system (sg b = -1) w has its asymptotic form for r >> 1 at any Q. *)
VerificationTest[
  QNMRadialKN[-1, 1, 1, 0, 0.3, 1.*^-10][{3., 10., 100.}]
  ,
  {1.5508988855985773624 + 1.4547517211808068695 I, 20.084120490370475015 - 16.382258524351073145 I,
   9.454330832936173958*^8 + 5.704934938860427439*^8 I}
  ,
  SameTest -> (Max[Abs[#1/#2 - 1]] < 10^-10 &),
  TestID -> "KNRadial-smallQ-EM"
]

(* For the gravitational-led family in the spin system -1 (sg b = -1); the frequency is passed from the default spin
   system, in which it is computed more accurately *)
VerificationTest[
  QNMRadialKN[-2, 2, 2, 0, 0.6, 1.*^-10, Method -> {"HintzSeparated", "SpinSystem" -> -1},
    "Frequency" -> QNMFrequencyKN[-2, 2, 2, 0, 0.6, 1.*^-10]][{3., 10., 100.}]
  ,
  {-0.82849673718302256216 + 0.92261708040247500215 I, 3.2325511169207129895 - 3.8591271636752165238 I,
   -1.6013088280142835448*^6 - 3.5649136526985546298*^7 I}
  ,
  SameTest -> (Max[Abs[#1/#2 - 1]] < 10^-10 &),
  TestID -> "KNRadial-smallQ-grav-spinsystem-1"
]

(* For sg b = 1 at small charge the expansion of w in 1/r holds only for r >> 1/Q^2; the function says so *)
VerificationTest[
  Head[QNMRadialKN[-2, 2, 2, 0, 0.6, 0.001]]
  ,
  QNMRadialFunction
  ,
  {QNMRadialKN::asym, QNMRadialKN::spinsys},
  TestID -> "KNRadial-smallQ-grav-message"
]


(* ::Section:: *)
(*Invalid input*)


VerificationTest[
  QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5, Method -> {"HintzSeparated", "Coordinates" -> "Kerr-Schild"}]
  ,
  $Failed
  ,
  {QNMRadialKN::coords},
  TestID -> "KNRadial-invalid-coordinates"
]

VerificationTest[
  QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5, Method -> {"HintzSeparated", "Component" -> "z"}]
  ,
  $Failed
  ,
  {QNMRadialKN::comp},
  TestID -> "KNRadial-invalid-component"
]

VerificationTest[
  QNMRadialKN[-2, 2, 2, 0, 0.6, 0]
  ,
  $Failed
  ,
  {QNMRadialKN::kerr},
  TestID -> "KNRadial-Q0"
]

VerificationTest[
  QNMRadialKN[0, 2, 2, 0, 0.6, 0.5]
  ,
  $Failed
  ,
  {QNMRadialKN::spin},
  TestID -> "KNRadial-scalar-not-covered"
]
