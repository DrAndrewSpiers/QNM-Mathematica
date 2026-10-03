(* ::Package:: *)

(* Tests for QNMFrequencyKN (Kerr-Newman quasinormal modes from Hintz's separated system).
   Run from a fresh kernel with
     TestReport["/path/to/QNM-Mathematica/Tests/QNMFrequencyKN.wlt"]
   or from a shell with
     wolframscript -code 'Print[TestReport["Tests/QNMFrequencyKN.wlt"]["TestsFailedCount"]]'
   Runtime: about a minute.

   Reference values:
   - RN and KN values: A. Spiers, "Quasinormal modes from Hintz's separated Kerr-Newman equations" (2026),
     Tables 1, 5, 6 (10 decimals; Chandrasekhar's RN equations agree with Table 1 to 2e-11).
   - Dias, Godazgar & Santos, PRL 114, 151101 (2015), data file Grav_KN_l2_a=Q.dat (50 digits; columns a = Q, Im, Re).
   - Mark, Yang, Zimmerman & Chen, PRD 91, 044025 (2015), O(Q^2) coefficient at a = 0.6, s = -2, l = m = 2.
   - RN overtone n = 1: 24-digit value computed with this package; the note quotes 0.355212988 - 0.275684377 I. *)

$QNMRepository = ParentDirectory[DirectoryName[$TestFileName]];
PacletDirectoryLoad[$QNMRepository];
Get["QNM`"];

VerificationTest[
  StringStartsQ[FindFile["QNM`"], $QNMRepository]
  ,
  True
  ,
  TestID -> "KN-LoadsRepositoryVersion"
]


(* ::Section:: *)
(*Reissner-Nordstrom (a = 0): Table 1 of the note*)


VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0., 0.5]
  ,
  0.3816771514 - 0.0896123793 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-RN-grav-l2-Q0.5"
]

VerificationTest[
  QNMFrequencyKN[-1, 1, 0, 0, 0., 0.5]
  ,
  0.2706902309 - 0.0950930956 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-RN-EM-l1-Q0.5"
]

VerificationTest[
  QNMFrequencyKN[-1, 2, 0, 0, 0., 0.8]
  ,
  0.5701302340 - 0.0990690623 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-RN-EM-l2-Q0.8"
]

VerificationTest[
  QNMFrequencyKN[-2, 3, 0, 0, 0., 0.3]
  ,
  0.6033253911 - 0.0929019229 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-RN-grav-l3-Q0.3"
]

(* small charge: the 3/(2Q) terms of Hintz's system must not cost accuracy *)
VerificationTest[
  QNMFrequencyKN[-2, 2, 0, 0, 0., 0.01]
  ,
  0.3736742664 - 0.0889625981 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-RN-grav-l2-Q0.01"
]

(* at a = 0 the frequency does not depend on m *)
VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0., 0.5] - QNMFrequencyKN[-2, 2, -1, 0, 0., 0.5]
  ,
  0.
  ,
  SameTest -> (Abs[#1 - #2] < 10^-10 &),
  TestID -> "KN-RN-m-independence"
]


(* ::Section:: *)
(*Kerr-Newman: Tables 5 and 6 of the note*)


VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.6, 0.5]
  ,
  0.5361271297 - 0.0811783720 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-grav-l2-m2-a0.6-Q0.5"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, -1, 0, 0.4, 0.7]
  ,
  0.3691169665 - 0.0891531928 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-grav-l2-m-1-a0.4-Q0.7"
]

VerificationTest[
  QNMFrequencyKN[-1, 2, 1, 0, 0.4, 0.7]
  ,
  0.5919044665 - 0.0939405660 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-EM-l2-m1-a0.4-Q0.7"
]

VerificationTest[
  QNMFrequencyKN[-1, 1, -1, 0, 0.6, 0.7]
  ,
  0.2575737779 - 0.0972079568 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-EM-l1-m-1-a0.6-Q0.7"
]


(* ::Section:: *)
(*Literature*)


(* Dias-Godazgar-Santos, rows 26 and 46 of Grav_KN_l2_a=Q.dat (a = Q) *)
VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, N[0.56534576476104024539], N[0.56534576476104024539]]
  ,
  0.53453334795195182650 - 0.08127868941517654718 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-10 &),
  TestID -> "KN-DGS-row26"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, N[0.70312591623312445921], N[0.70312591623312445921]]
  ,
  0.83334200859157986864 - 0.02929026034870258779 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-10 &),
  TestID -> "KN-DGS-row46"
]

(* ::Section:: *)
(*Near extremality*)


(* Reissner-Nordstrom, Q = 0.99: Chandrasekhar benchmark quoted in Section 8 of the note *)
VerificationTest[
  QNMFrequencyKN[-2, 2, 0, 0, 0., 0.99]
  ,
  0.429296793341 - 0.084265944127 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-10 &),
  TestID -> "KN-RN-near-extremal-Q0.99"
]

(* a = Q at 1 - a^2 - Q^2 = 1.5*10^-5 (row 59 of the Dias-Godazgar-Santos file, whose value differs by 7*10^-8);
   reference computed with this package at 24 digits *)
VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, N[0.70710161483605208157], N[0.70710161483605208157]]
  ,
  0.93840157711047367393 - 0.00126540082026782121 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-10 &),
  TestID -> "KN-near-extremal-a=Q"
]


(* ::Section:: *)
(*Small charge*)


(* Mark-Yang-Zimmerman-Chen: omega = omega_Kerr + Q^2 omega1 + O(Q^4); omega1 at a = 0.6 for s = -2, l = m = 2 *)
VerificationTest[
  Module[{a = 0.6, qs = {0.2, 0.1, 0.05}, wK = 0.49404478178138433 - 0.08376520216104169 I, oms},
    oms = QNMFrequencyKN[-2, 2, 2, 0, a, #] & /@ qs;
    First[LinearSolve[Table[{1, q^2, q^4}, {q, qs}], (oms - wK)/qs^2]]]
  ,
  0.12895437458774348 + 0.004480822088666132 I
  ,
  SameTest -> (Abs[#1 - #2] < 5 10^-7 &),
  TestID -> "KN-MYZC-omega1-a0.6"
]


(* ::Section:: *)
(*Overtones*)


VerificationTest[
  QNMFrequencyKN[-2, 2, 0, 1, 0., 0.5]
  ,
  0.355212987535554910936 - 0.275684376707070168995 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-RN-grav-l2-n1-Q0.5"
]

(* the KN overtone tends to the Kerr overtone as Q -> 0 (the difference is O(Q^2));
   the Kerr reference may emit FindRoot::lstol, which is not part of this test *)
VerificationTest[
  Abs[QNMFrequencyKN[-2, 2, 2, 1, 0.5, 0.01] - Quiet[QNMFrequency[-2, 2, 2, 1, 0.5]]] < 10^-4
  ,
  True
  ,
  TestID -> "KN-Kerr-limit-n1"
]


(* ::Section:: *)
(*Limits, symmetries and precision*)


VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.3, 0]
  ,
  QNMFrequency[-2, 2, 2, 0, 0.3]
  ,
  TestID -> "KN-Q0-is-Kerr"
]

VerificationTest[
  {QNMFrequencyKN[-2, 2, 2, 0, 0.3, -0.5], QNMFrequencyKN[2, 2, 2, 0, 0.3, 0.5], QNMFrequencyKN[-2, 2, -2, 0, -0.3, 0.5]}
  ,
  ConstantArray[QNMFrequencyKN[-2, 2, 2, 0, 0.3, 0.5], 3]
  ,
  SameTest -> (Max[Abs[#1 - #2]] < 10^-10 &),
  TestID -> "KN-symmetries"
]

VerificationTest[
  Module[{w = QNMFrequencyKN[-2, 2, 0, 0, 0, 0.5`24]},
    {Precision[w] >= 23, Abs[w - (0.381677151370047620173437560063 - 0.089612379287727120475961600335 I)] < 10^-22}]
  ,
  {True, True}
  ,
  TestID -> "KN-arbitrary-precision"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, {1, 2}, 0, 0.3, 0.5] // Length
  ,
  2
  ,
  TestID -> "KN-Listable"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.6, 0.5, Method -> {"HintzSeparated", "SeriesOrder" -> 30, "MatchingRadius" -> 5, "RayLength" -> 60}]
  ,
  0.5361271297 - 0.0811783720 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-method-options"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.3, 0.5, Method -> {"HintzSeparated", "InitialGuess" -> 0.436 - 0.088 I}]
  ,
  0.43633377254739938830 - 0.08783936104152728126 I
  ,
  SameTest -> (Abs[#1 - #2] < 10^-9 &),
  TestID -> "KN-initial-guess"
]


(* ::Section:: *)
(*Invalid input*)


VerificationTest[
  QNMFrequencyKN[-2, 1, 0, 0, 0.3, 0.5]
  ,
  $Failed
  ,
  {QNMFrequencyKN::params},
  TestID -> "KN-invalid-l"
]

VerificationTest[
  QNMFrequencyKN[0, 2, 0, 0, 0.3, 0.5]
  ,
  $Failed
  ,
  {QNMFrequencyKN::spin},
  TestID -> "KN-scalar-not-covered"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.8, 0.6]
  ,
  $Failed
  ,
  {QNMFrequencyKN::extremal},
  TestID -> "KN-extremal"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.3 + 0.1 I, 0.5]
  ,
  $Failed
  ,
  {QNMFrequencyKN::cmplx},
  TestID -> "KN-complex-spin"
]

VerificationTest[
  QNMFrequencyKN[-2, 2, 2, 0, 0.3, 0.5, Method -> "NoSuchMethod"]
  ,
  $Failed
  ,
  {QNMFrequencyKN::optx},
  TestID -> "KN-unknown-method"
]
