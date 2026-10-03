(* ::Package:: *)

(* ::Title:: *)
(*QNM*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["QNM`",
  {
  "SpinWeightedSpheroidalHarmonics`",
  "Teukolsky`",
  "Teukolsky`TeukolskyRadial`"
  }
];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{QNMFrequency, QNMFrequencyKN, QNMRadial, QNMRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


QNMFrequency::usage = "QNMFrequency[s, l, m, n, a] computes a quasinormal mode frequency.";


QNMFrequencyKN::usage = "QNMFrequencyKN[s, l, m, n, a, Q] computes a quasinormal mode frequency of a Kerr-Newman black hole with spin a and charge Q (in units of the mass M). \
The mode is the coupled gravito-electromagnetic mode that reduces to the Kerr mode QNMFrequency[s, l, m, n, a] as Q \[Rule] 0: \
|s| = 2 selects the gravitational-led and |s| = 1 the electromagnetic-led family. It is computed from P. Hintz's separation of the Einstein-Maxwell perturbation equations.";


QNMRadial::usage = "QNMRadial[s, l, m, n, a] computes the radial eigenfunction of a quasinormal mode.";


QNMRadialFunction::usage = "QNMRadialFunction[...] is an object representing a quasinormal mode solution to the radial Teukolsky equation.";


(* ::Subsection::Closed:: *)
(*Error messages*)


QNMFrequency::nointerp = "Interpolation data not available for s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::optx = "Unknown options in `1`.";
QNMFrequency::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequency::findroot = "FindRoot failed to converge to the requested accuracy.";
QNMFrequency::cmplx = "Only real values of a are allowed, but a=`1` specified.";
QNMFrequency::nokerr = "Method \"`1`\" only supported for Schwarzschild spacetime, but a=`2` specified.";
QNMFrequency::acc = "Accuracy of the calculated quasinormal mode frequency `1` is lower than that of the initial guess `2`.";
QNMFrequency::noacc = "Failed estimating accuracy of the calculated quasinormal mode frequency `1` compared to the initial guess `2`.";
QNMFrequency::allln = "Method `1` can only be used with l=All and n=All to compute a set of frequencies but l=`2` and n=`3` specified.";
QNMFrequencyKN::optx = "Unknown options in `1`.";
QNMFrequencyKN::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMFrequencyKN::spin = "Kerr-Newman quasinormal modes are only available for the coupled gravito-electromagnetic perturbations, |s| = 1 (electromagnetic-led) or |s| = 2 (gravitational-led), but s=`1` specified.";
QNMFrequencyKN::cmplx = "Only real values of a and Q are allowed, but a=`1` and Q=`2` specified.";
QNMFrequencyKN::extremal = "The spin a=`1` and charge Q=`2` must satisfy a^2 + Q^2 < 1 (sub-extremal black hole).";
QNMFrequencyKN::seed = "Failed to compute the Kerr frequency QNMFrequency[`1`, `2`, `3`, `4`, `5`] used to start the continuation in Q.";
QNMFrequencyKN::cont = "Continuation of the quasinormal mode from the Kerr value failed at a=`1`, Q=`2`.";
QNMFrequencyKN::findroot = "The root search starting from `1` did not converge.";
QNMFrequencyKN::acc = "The frequency `1` is estimated to be accurate only to about `2` (working precision `3`). Near extremality increase the precision of a and Q or set Method -> {\"HintzSeparated\", WorkingPrecision -> 30}.";
QNMRadial::optx = "Unknown options in `1`.";
QNMRadial::params = "Invalid parameters s=`1`, l=`2`, m=`3`, n=`4`.";
QNMRadial::coords = "Coordinate options are either \"BL\", \"Boyer-Lindquist\", or \"Hyperboloidal\", but got `1`.";
QNMRadial::convergence = "Eigenvalue failed to converge to specified tolerance. Final value `1`.";
QNMRadialFunction::dmval = "Radius `1` lies outside the computational domain.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Global Debug Flag*)


DEBUG=False;


(* ::Section::Closed:: *)
(*Utilities*)


(* ::Subsection::Closed:: *)
(*Horizon locations*)


M=1;


rp[a_, M_] := M+Sqrt[M^2-a^2];
rm[a_, M_] := M-Sqrt[M^2-a^2];


(* ::Subsection::Closed:: *)
(*Discretised Teukolsky operator on a hyperboloidal slice*)


\[ScriptCapitalM][s_, m_, a_, \[Omega]_, \[Lambda]_, n_] :=
 Module[{\[Rho], \[CapitalDelta], R, dR, d2R, A, B, M=1},
  \[Rho] = 1/rp[a, M] Reverse[1/2 (1+Cos[\[Pi] Subdivide[n-1]])];
  \[CapitalDelta] = 1-2M \[Rho]+a^2 \[Rho]^2;
  R = IdentityMatrix[n];
  dR  = NDSolve`FiniteDifferenceDerivative[Derivative[1], \[Rho], DifferenceOrder -> "Pseudospectral", PeriodicInterpolation -> False]["DifferentiationMatrix"];
  d2R = NDSolve`FiniteDifferenceDerivative[Derivative[2], \[Rho], DifferenceOrder -> "Pseudospectral", PeriodicInterpolation -> False]["DifferentiationMatrix"];
  A = 2 I \[Omega]-2(1+s)\[Rho]+2(I \[Omega](a^2-8M^2)+I m a +(s+3)M)\[Rho]^2+4(2 I \[Omega] M-1) a^2 \[Rho]^3;
  B = (a^2-16M^2)\[Omega]^2+2(m a +2 I s M)\[Omega]+2(4(a^2-4M^2)M \[Omega]^2+(4 m a M-4I (s+2) M^2+I a^2)\[Omega]+ I m a+(s+1)M)\[Rho]+2(8 M^2 \[Omega]^2+6 I M \[Omega]-1)a^2 \[Rho]^2;
  -\[Rho]^2 \[CapitalDelta] d2R + A dR + (B+(\[Lambda] + 2 a m \[Omega] - (a \[Omega])^2)) R
]


(* ::Subsection::Closed:: *)
(*Discretised 2D Teukolsky operator on a hyperboloidal slice*)


\[ScriptCapitalM]2D[s_, m_, a_, nx_, nz_] :=
 Module[{\[Chi], z0, z1, \[CapitalDelta]z, z, x, Id, Zero, Dx, D2x, D\[Sigma], D2\[Sigma], L1\[Sigma]\[Sigma], L1\[Sigma], L10, L1xx, L1x, L2\[Sigma], L20, W, L1, L2,
         rh, \[Kappa], \[Stigma], \[Lambda], \[Delta]1, \[Delta]2, \[Sigma], X},
  (* Constants *)
  rh = rp[a, M];
  \[Kappa] = a/rh;
  \[Stigma] = s;
  \[Lambda] = M;
  \[Delta]1 = Abs[m-s];
  \[Delta]2 = Abs[m+s];

  (* Radial and angular coordinates *)
  z0 = 0;
  z1 = 1;
  \[CapitalDelta]z = z1-z0;
  \[Chi] = N@Reverse[Cos[\[Pi] Subdivide[nz-1]]];
  z = z0+1/2 \[CapitalDelta]z*(1+\[Chi]);
  x = N@Reverse[Cos[\[Pi] Subdivide[nx-1]]];

  (* Discretised operators *)
  Id = IdentityMatrix[nx*nz];
  Zero = ConstantArray[0, {nx*nz, nx*nz}];
  Dx  = NDSolve`FiniteDifferenceDerivative[Derivative[0,1], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];
  D2x = NDSolve`FiniteDifferenceDerivative[Derivative[0,2], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];
  D\[Sigma]  = NDSolve`FiniteDifferenceDerivative[Derivative[1,0], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];
  D2\[Sigma] = NDSolve`FiniteDifferenceDerivative[Derivative[2,0], {z, x}, DifferenceOrder -> {"Pseudospectral", "Pseudospectral"}, PeriodicInterpolation -> {False, False}]["DifferentiationMatrix"];

  (* Coordinates on the full (flattened) 2D grid *)
  {\[Sigma], X} = {Flatten[Outer[#1&, z, x]], Flatten[Outer[#2&, z, x]]};

  L1xx = 1-X^2;
  L1x = \[Delta]1-\[Delta]2-X (2+\[Delta]1+\[Delta]2);
  
  L1\[Sigma]\[Sigma] = (-1+\[Sigma]) \[Sigma]^2 (-1+\[Kappa]^2 \[Sigma]);
  L1\[Sigma] = \[Sigma] (4 \[Kappa]^2 \[Sigma]^2+2 (1+\[Stigma])-\[Sigma] (2 I m \[Kappa]+(1+\[Kappa]^2) (3+\[Stigma])));
  L10 = 1/2 (-m^2-\[Delta]2-\[Delta]1 (1+\[Delta]2)-4 I m \[Kappa] \[Sigma]+4 \[Kappa]^2 \[Sigma]^2-2 (1+\[Kappa]^2) \[Sigma] (1+\[Stigma])+\[Stigma] (2+\[Stigma]));
  
  L2\[Sigma] = (2 rh (1+\[Sigma]^2 (-2+2 \[Kappa]^4 (-1+\[Sigma])+\[Kappa]^2 (-3+2 \[Sigma]))))/\[Lambda];
  L20 = (2 rh (3 (\[Kappa]^2+\[Kappa]^4) \[Sigma]^2-I m (\[Kappa]+2 \[Kappa] (1+\[Kappa]^2) \[Sigma])+\[Stigma]+\[Kappa] (-I X+\[Kappa]) \[Stigma]-\[Sigma] (2+\[Stigma]+\[Kappa]^4 (2+\[Stigma])+\[Kappa]^2 (3+2 \[Stigma]))))/\[Lambda];
  
  W = (rh^2 (4 (1+\[Sigma])+\[Kappa]^2 (7+X^2+4 \[Kappa]^2+4 (2+2 \[Kappa]^2+\[Kappa]^4) \[Sigma]-4 (1+\[Kappa]^2)^2 \[Sigma]^2)))/\[Lambda]^2;
  
  L1 = L1\[Sigma]\[Sigma] D2\[Sigma] + L1\[Sigma] D\[Sigma] + L10 Id + L1xx D2x + L1x Dx;
  L2 = L2\[Sigma] D\[Sigma] + L20 Id;
  ArrayFlatten[{{Zero, Id}, {L1/W, L2/W}}]
];


(* ::Subsection::Closed:: *)
(*Chebyshev interpolation*)


(* Convert from values at Chebyshev nodes to Chebyshev coefficients *)
chebCoeffs[values_] :=
 Module[{coeffs},
  coeffs = Sqrt[2/(Length[values]-1)] FourierDCT[values, 1];
  coeffs[[{1,-1}]] /= 2;
  coeffs
];


(* Construct a Chebyshev InterpolatingFunction *)
(* Useful references on the structure of an InterpolatingFunction:
   https://mathematica.stackexchange.com/a/28341
   https://mathematica.stackexchange.com/a/98349
   https://mathematica.stackexchange.com/a/114065 *)
chebInterp[data_, domain_] :=
 Module[{order, coeffs},
  order = Length[data] - 1;
  coeffs = chebCoeffs[data];
  InterpolatingFunction[
    {domain}    (* Interpolation domain *),
    {5          (* InterpolatingFunction version *),
     1,         (* Bit field: always 1 for Chebyshev *)
     order,     (* Max derivative order in data *)
     {2},       (* Domain grid size *)
     {order+1}, (* Interpolation order + 1 *)
     0,         (* Not a derivative of an existing InterpolatingFunction *)
     0,         (* Non-periodic interpolation *)
     0, 0,
     Automatic, (* Extrapolation handler *)
     {}, {}, False},
    {domain},   (* Domain grid *)
    {data, coeffs},
    {{{{1}, {1,2}, {1,2}}, {Automatic, ChebyshevT, ChebyshevT}}}
  ]
];


(* ::Subsection::Closed:: *)
(*Functions for Leaver's method*)


\[Beta][s_] := s^2 -1;
k1[m_, s_] := 1/2 Abs[m-s];
k2[m_, s_] := 1/2 Abs[m+s];


(* ::Subsubsection::Closed:: *)
(*Schwarzschild*)


(* Functions for the Continued Fraction used in Leaver's method for Schwarzschild *)
\[Alpha][i_, \[Omega]_]  := i^2 - 4 M I \[Omega] i + 2 i - 4 M I \[Omega] + 1;
\[Delta][i_, \[Omega]_,s_, l_] := -2 i^2 - (2 - 16 M I \[Omega])i + 32 M^2 \[Omega]^2 + 8 M I \[Omega] - l(l + 1) + \[Beta][s];
\[Gamma][i_, \[Omega]_, s_] := i^2 - 8 M I \[Omega] i - 16 M^2 \[Omega]^2 - \[Beta][s] - 1;


(* ::Subsubsection::Closed:: *)
(*Kerr*)


b[a_] := Sqrt[4 M^2 - 4 a^2];


(* Functions for the Continued Fraction used in Leaver's method for Kerr *)
\[Alpha]freq[i_,\[Omega]_, s_, m_, a_]:= i^2 + (2-s-2 M I \[Omega] -2 I/b[a] (2 M^2 \[Omega] - a m))i + 1 - s - 2 M I \[Omega] -2 I/b[a] (2 M^2 \[Omega] - a m);
\[Beta]freq[i_,\[Omega]_, Alm_, s_, m_, a_]:= -2 i^2 +(-2 + 2(4 M + b[a]) I \[Omega] + 4 I/b[a] (2 M^2 \[Omega] - a m)) i + (16 M^2 + 4 M b[a] - a^2) \[Omega]^2 - s - 1 - 2 a m \[Omega] - Alm +(4 M + b[a]) I \[Omega] + (8 M \[Omega] + 2 I)/b[a] (2 M^2 \[Omega] - a m);
\[Gamma]freq[i_, \[Omega]_, s_, m_, a_]:= i^2 + (s - 6 M I \[Omega] - 2 I/b[a] (2 M^2 \[Omega] - a m)) i - 8 M^2 \[Omega]^2 - 4 M I \[Omega] s - 8 M \[Omega]/b[a] (2 M^2 \[Omega] - a m);

\[Alpha]ang[i_, \[Omega]_, s_, m_, a_]:= -2 (i + 1) (i + 2 k1[m, s] + 1);
\[Beta]ang[i_, \[Omega]_, Alm_, s_, m_, a_]:= i (i - 1) + 2 i (k1[m, s] + k2[m, s] + 1 - 2 a \[Omega]) - (2 a \[Omega] (2 k1[m, s] + s + 1) - (k1[m, s] + k2[m, s]) (k1[m, s] + k2[m, s] + 1)) - (a^2 \[Omega]^2 + s(s+1) + Alm);
\[Gamma]ang[i_, \[Omega]_, s_, m_, a_]:= 2 a \[Omega] (i + k1[m, s] + k2[m, s] + s);


(* ::Subsection::Closed:: *)
(*Continued fraction*)


(* ::Text:: *)
(*This is preferred over Mathematica's ContinedFractionK function as we can get an error estimate on the result using this function.*)


CF[a_, b_, {n_, n0_}] := 
  Module[{A, B, ak, bk, res = Indeterminate, j = n0},
   A[n0 - 2] = 1;
   B[n0 - 2] = 0;
   ak[k_] := ak[k] = (a /. n -> k);
   bk[k_] := bk[k] = (b /. n -> k);
   A[n0 - 1] = 0(*bk[n0-1]*);
   B[n0 - 1] = 1;
   A[k_] := A[k] = bk[k] A[k - 1] + ak[k] A[k - 2];
   B[k_] := B[k] = bk[k] B[k - 1] + ak[k] B[k - 2];
   While[Quiet[Check[Abs[1-(A[j-1]/B[j-1])/(A[j]/B[j])], j--; False, General::munfl], General::munfl] > 10^(1-Precision[A[j]/B[j]]), j++];
   res = A[j]/B[j];
   Clear[A, B, ak, bk];
   res
];


(* ::Section::Closed:: *)
(*QNMFrequency*)


prec[a_] :=
 Which[
   a != 0,
   Precision[a],
   a == 0 && Precision[a] == MachinePrecision,
   MachinePrecision,
   a == 0,
   Accuracy[a]
 ];


(* ::Subsection::Closed:: *)
(*Asymptotic approximations*)


(* ::Text:: *)
(*Multiple asymptotic expansions are used for Schwarzschild, with different expansions providing better approximations for different values of l and n.*)


(* ::Subsubsection::Closed:: *)
(*Schwarzschild asymptotic expansion for l >> n, l>> 1*)


(* ::Text:: *)
(*This expansion is taken from Dolan, Ottewill, Classical and Quantum Gravity, Vol. 26, 2009.*)


Schwarzfinit1[s_, l_, n_]:= ((I*(n+1/2)*(\[Beta][s]^2/27 + (\[Beta][s]*(1100*(n+1/2)^2 - 2719))/46656 + (11273136*(n+1/2)^4 - 52753800*(n+1/2)^2 + 66480535)/2902376448))/(l+1/2)^4 + 
   (-(\[Beta][s]^2/27) + (\[Beta][s]*(204*(n+1/2)^2 + 211))/3888 + (854160*(n+1/2)^4 - 1664760*(n+1/2)^2 - 776939)/40310784)/(l+1/2)^3 - (I*(n+1/2)*(\[Beta][s]/9 + (235*(n+1/2)^2)/3888 - 1415/15552))/(l+1/2)^2 + 
   (\[Beta][s]/3 - (5*(n+1/2)^2)/36 - 115/432)/(l+1/2) + (l+1/2) - I*(n+1/2))/(Sqrt[27]*M);


(* ::Subsubsection::Closed:: *)
(*Schwarzschild low order asymptotic expansion for n>>l, n>>1*)


(* ::Text:: *)
(*The expansion is taken from Casals, Dolan, Ottewill, Wardell, Phys. Rev. D, Vol. 88, 2013*)


Schwarzfinit2[n_] := Log[3]/(8 \[Pi] M) - I (n+1/2)/(4 M);


(* ::Subsubsection::Closed:: *)
(*Kerr*)


(* ::Text:: *)
(*For Kerr, an initial guess is needed for both the frequency and the spheroidal eigenvalue Subscript[A, lm].*)


Kerrfinit[s_, l_, m_, n_, a_] :=
 Module[{b, \[CapitalDelta], \[Mu], Eikonal, Rp, \[CapitalOmega]r, \[CapitalOmega]i, finit},
  \[CapitalDelta][r_] := r^2 -2 M r + a^2;
  \[Mu] = m/(l+1/2);(* Useful parameter *)
  Eikonal[rp_]:= 2(rp/M)^4(rp/M - 3)^2 + 4 (rp/M)^2((1 - \[Mu]^2)(rp/M)^2 - 2(rp/M) - 3(1 - \[Mu]^2))(a/M)^2 + (1 - \[Mu]^2) ( (2 - \[Mu]^2) (rp/M)^2 + 2 (2 + \[Mu]^2) (rp/M) + (2 - \[Mu]^2)) (a/M)^4;

  Rp = FindRoot[Eikonal[rp]==0, {rp, 3}] [[1]][[2]]; (* FIXME: This always produces a machine-precision result. Maybe use Solve instead? *)

  \[CapitalOmega]r = -If[\[Mu] == 0, \[Pi]/2 Sqrt[\[CapitalDelta][Rp]]/((Rp^2 + a^2) EllipticE[(a^2 \[CapitalDelta][Rp])/((Rp^2 + a^2)^2)] ) , (M - Rp) \[Mu] a / ((Rp - 3M)Rp^2 + (Rp + M) a^2)];
  \[CapitalOmega]i = \[CapitalDelta][Rp](Sqrt[4(6Rp^2 \[CapitalOmega]r^2 -1) + 2a^2\[CapitalOmega]r^2(3 - \[Mu]^2)])/(2 Rp^4 \[CapitalOmega]r - 4 a M Rp \[Mu] + a^2 Rp \[CapitalOmega]r(Rp(3 - \[Mu]^2) + 2 M(1+\[Mu]^2)) + a^4\[CapitalOmega]r(1-\[Mu]^2));

  finit = (l + 1/2) Abs[\[CapitalOmega]r] - I (n + 1/2) Abs[\[CapitalOmega]i]
];


KerrAinit[s_, l_, m_, n_, a_] := (l+1/2)^2 - (a Kerrfinit[s, l, m, n, a])^2 /2 (1 - (m/(l+1/2))^2);


(* ::Subsection::Closed:: *)
(*Interpolation of tabulated QNM frequencies*)


$QNMInstallationDirectory = FileNameDrop[FindFile["QNM`"], -2];
$QNMDataDirectory = FileNameJoin[{$QNMInstallationDirectory, "Data"}];


ClearAll[QNMFrequencyInterpolation];

Options[QNMFrequencyInterpolation] = Options[Interpolation];

QNMData[file_, dataset_] := QNMData[file, dataset] = Quiet[Import[file, {"Datasets", dataset}, "ComplexKeys"->{"r", "i"}], {Import::general, Import::noelem, Import::nffil}];

QNMFrequencyInterpolation[s_, l_, m_, n_, opts:OptionsPattern[]] :=
 Module[{h5file, dataset, data, ret},
  h5file = FileNameJoin[{$QNMDataDirectory, "QNM_s"<>ToString[s]<>".h5"}];
  dataset = "/l"<>ToString[l]<>"/m"<>ToString[m]<>"/n"<>ToString[n];
  data = QNMData[h5file, dataset];
  If[MatchQ[data, <|"a"->_, "omega"->_|>],
    ret = Interpolation[Transpose[Lookup[data, {"a","omega"}]], opts];,
    Message[QNMFrequency::nointerp, s, l, m, n];
    ret = Function[{a}, $Failed];
  ];
  ret
];


(* ::Subsection::Closed:: *)
(*Calculate QNM frequency by finding zero of incidence amplitude*)


Options[QNMFrequencyInIncidenceAmplitude] = {"InitialGuess" -> Automatic};


QNMFrequencyInIncidenceAmplitude[s_Integer, l_Integer, m_Integer, n_, a_, OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, inInc, prec = prec[a]},
  \[Omega]guess=OptionValue["InitialGuess"];
  If[\[Omega]guess === Automatic,
    \[Omega]guess = Quiet[QNMFrequencyInterpolation[s, l, m, n][a], QNMFrequency::nointerp];
  ];
  If[\[Omega]guess == $Failed, 
    Message[QNMFrequency::nointerp, s, l, m, n];
    If[a==0,
      Which[
        l <= Max[n,2], 
        \[Omega]guess = SpectralInitialGuess[s, l, n];,
        l>n,
        \[Omega]guess = Schwarzfinit1[s, l, n];
      ];,
      \[Omega]guess = Kerrfinit[s, l, m, n, a];
    ];
  ];
  inInc[\[Omega]_?NumericQ] := inInc[\[Omega]] = TeukolskyRadial[s, l, m, If[a==0, 0, a], \[Omega], Method -> "MST"]["In"]["Amplitudes"]["Incidence"];
  \[Omega]QNM /. Quiet[Check[FindRoot[inInc[\[Omega]QNM], {\[Omega]QNM, SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec], \[Omega]QNM -> $Failed, FindRoot::nlnum], FindRoot::nlnum]
];


(* ::Subsection::Closed:: *)
(*Ripley Hyperboloidal method*)


Options[QNMFrequencyHyperboloidal] = {"NumPoints" -> 32, "InitialGuess" -> Automatic, "AccuracyCheck" -> True};


QNMFrequencyHyperboloidal[s_Integer, l_Integer, m_Integer, n_Integer, a_, opts:OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, \[Omega], \[Delta]\[Lambda], numpoints, prec = prec[a]},
  numpoints = OptionValue["NumPoints"];
  \[Omega]guess=OptionValue["InitialGuess"];
  If[\[Omega]guess === Automatic,
    \[Omega]guess = Quiet[QNMFrequencyInterpolation[s, l, m, n][a], QNMFrequency::nointerp];
  ];
  If[\[Omega]guess == $Failed, 
    Message[QNMFrequency::nointerp, s, l, m, n];
    If[a==0,
      Which[
        l <= Max[n,2], 
        \[Omega]guess = SpectralInitialGuess[s, l, n];,
        l>n,
        \[Omega]guess = Schwarzfinit1[s, l, n];
      ];,
      \[Omega]guess = Kerrfinit[s, l, m, n, a];
    ];
  ];
  \[Delta]\[Lambda][\[Omega]_?NumericQ] := \[Delta]\[Lambda][\[Omega]] = Module[{\[Lambda], Mat},
    \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
    Mat = \[ScriptCapitalM][s, m, a, \[Omega], \[Lambda], numpoints];
    Eigenvalues[Mat, -1]
  ];
  \[Omega]QNM = \[Omega] /. FindRoot[\[Delta]\[Lambda][\[Omega]], {\[Omega], SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec];

  If[OptionValue["AccuracyCheck"] == True &&
     Quiet[Check[Abs[TeukolskyRadial[s, l, m, If[a==0, 0, a], \[Omega]QNM, Method -> "MST"]["In"]["Amplitudes"]["Incidence"] / 
               TeukolskyRadial[s, l, m, If[a==0, 0, a], \[Omega]guess, Method -> "MST"]["In"]["Amplitudes"]["Incidence"]] > 10^3,
           Message[QNMFrequency::noacc, \[Omega]QNM, \[Omega]guess]; False], All, QNMFrequency::noacc],
    Message[QNMFrequency::acc, \[Omega]QNM, \[Omega]guess];
  ];

  \[Omega]QNM
];


(* ::Subsection::Closed:: *)
(*Ansorg-Macedo hyperboloidal method*)


(* ::Text:: *)
(*Spectral method using hyperboloidal slicing. This technique was adapted from a Reissner-Nordstrom version kindly provided by Rodrigo Macedo, Queen Mary University of London. The technique is outlined in detail in Ansorg, Macedo, Phys. Rev. D, Vol. 93, 2016*)


SpectralInitialGuess[s_, l_, n_] :=
 Module[{L1a, L1b, L2a, L2b, L2c, Prec, Ndiv, ndiv, \[Sigma], \[Sigma]0, \[Sigma]1, \[CapitalDelta]\[Sigma], x, Id, Z, c, d\[Sigma], D\[Sigma], D2\[Sigma], L1, L2, M, Eigens, Filtered, sPattern},
  (* Operators appearing in the wave equation in coordinates (\[Tau], \[Sigma]), excluding radial derivatives *)
  L1a = -((2 \[Sigma])/(\[Sigma]+1));
  L1b = (1-2\[Sigma]^2)/(\[Sigma]+1);
  L2a = (-l(l+1) - \[Sigma](1-s^2))/(\[Sigma]+1);
  L2b = (\[Sigma](2-3\[Sigma]))/(\[Sigma]+1);
  L2c = (\[Sigma]^2 (1-\[Sigma]))/(\[Sigma]+1);

  Prec = 100; (*Numerical precision, machine precision found to be insufficient*)
	  		(* Make optional for users to set? Requires more testing as well *)
  If[n > 2, Ndiv = 6n, Ndiv = 20]; (* Subdivision of radial domain \[Sigma] \[Element] [0,1] used for discretization of radial derivs *)
  ndiv = Ndiv + 1; (* Matrix dimension *)

  \[Sigma]0 = 0;
  \[Sigma]1 = 1;
  \[CapitalDelta]\[Sigma] = \[Sigma]1-\[Sigma]0;

  (* The method for discretizing the radial derivatives is detailed in Spectral Methods in Matlab, Lloyd N. Trefethen *)
  x[i_]:= Cos[(i \[Pi])/Ndiv]; (* Chebyshev points *)
  \[Sigma] = N[Table[\[Sigma]0 + \[CapitalDelta]\[Sigma]  (1+x[i])/2, {i, 0, Ndiv}], Prec]; (* Radial coordinate grid *)

  Id = IdentityMatrix[ndiv];
  Z = ConstantArray[0, {ndiv, ndiv}]; (* Zero matrix *)
  c[i_] := If[i(i-Ndiv) ==0, 2, 1];

  (* Radial derivative using Chebyshev-Lobatto method *)
  d\[Sigma][i_, j_] := 2/\[CapitalDelta]\[Sigma] Which[i==j && i == 0, (2 Ndiv^2+1)/6, i==j && i == Ndiv, -((2 Ndiv^2+1)/6), i==j,  -x[i]/(2(1-x[i]^2)) , i!=j, c[i]/c[j] (-1)^(i+j)/(x[i]-x[j])] ;

  (* First deriv matrix *)
  D\[Sigma]=N[Table[d\[Sigma][i, j], {i, 0, Ndiv}, {j, 0, Ndiv}], Prec];
  (* Second deriv matrix *)
  D2\[Sigma] = D\[Sigma] . D\[Sigma];

  (* Full opertors in wave eq *)
  L1 = L1b D\[Sigma] + L1a Id;
  L2 = L2c D2\[Sigma] + L2b D\[Sigma] + L2a Id;

  M = ArrayFlatten[{{Z, Id}, {L2, L1}}];

  (* Solve for eigenvalues *)
  Eigens = Eigenvalues[M] ;

  sPattern = Which[s==0,x_/;(Im[x]==0. ), 
    Abs[s]==1, x_/;(Im[x]==0. ) ,
    Abs[s]==2,  x_/;(Im[x]==0. && Abs[Re[x]+4] >0.2 )]; 

  (* Remove all values with Im[\[Omega]] = 0, which actually correspond to the purely imaginary roots (see multiplication by I/4 below).
     Also remove those with Im[\[Omega]] > 0, as these are the corresponding QNMs in the 3rd quadrant. *)
  Filtered = Reverse[DeleteCases[Eigens,x_/;(Im[x]>=0. )]];

  (*Pick out eigenvalue corresponding to the desired overtone *)
  (*Would be nice to save this and use the other elements for initial seeds for other values of n *)
  (*However, the accuracy decreases as one goes down the list *)
  Which[Abs[s]==2 && l==2 && n>8, I/4 Filtered[[n]], Abs[s]==2 && l==2 && n==8, (0.4615178773933189 10^-15 - I 3.9999999999996)/2, True, I/4 Filtered[[n+1]]]
];


(* ::Subsection::Closed:: *)
(*Assaad-Macedo hyperboloidal method in Kerr*)


(* ::Text:: *)
(*Spectral method using hyperboloidal slicing in Kerr. This technique was adapted from a version kindly provided by Rodrigo Macedo. The technique is outlined in detail in Assaad and Macedo, arXiv:2506.04326.*)


Options[SpectralInitialGuessKerr] = {"NumAngularPoints" -> 10, "NumRadialPoints" -> 10};


SpectralInitialGuessKerr[s_, m_, a_, opts:OptionsPattern[]] :=
 Module[{nAng, nRad},
  nAng = OptionValue["NumAngularPoints"];
  nRad = OptionValue["NumRadialPoints"];
  
  ReverseSortBy[I Eigenvalues[\[ScriptCapitalM]2D[s, m, a, nAng, nRad]], Im]
];


(* ::Subsection::Closed:: *)
(*Leaver method*)


(* ::Text:: *)
(*Calculation based on the method in Leaver, Proc. R. Soc. Lond. A. 402, 1985, Nollert, Phys. Rev. D, Vol. 47, 1993 as well as the code provided online by Emanuele Berti, https://pages.jh.edu/~eberti2/ringdown/.*)


(* Equations of which the QNMs and spheroidal eigenvalues are roots *)
Leaver[\[Omega]_?NumericQ, s_?IntegerQ, l_?IntegerQ, nInv_?IntegerQ] :=
 Module[{n},
  \[Delta][nInv, \[Omega], s, l] + ContinuedFractionK[-\[Alpha][nInv-n, \[Omega]] \[Gamma][nInv-n+1, \[Omega], s], \[Delta][nInv-n,\[Omega],s, l], {n,1,nInv}] + CF[-\[Alpha][n-1, \[Omega]] \[Gamma][n, \[Omega], s], \[Delta][n, \[Omega], s, l], {n, nInv+1}]
];
  
Leaver31[\[Omega]_?NumericQ, s_?IntegerQ, l_?IntegerQ, m_?IntegerQ, a_?NumericQ, nInv_?IntegerQ] := Module[{Alm, n},
  Alm = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]] + 2 a m \[Omega] - a^2 \[Omega]^2;
  \[Beta]freq[nInv, \[Omega], Alm, s, m, a] + ContinuedFractionK[-\[Alpha]freq[nInv-n, \[Omega], s, m, a] \[Gamma]freq[nInv-n+1, \[Omega], s, m, a], \[Beta]freq[nInv-n,\[Omega], Alm, s, m, a], {n, 1, nInv}] + CF[-\[Alpha]freq[n-1, \[Omega], s, m, a] \[Gamma]freq[n, \[Omega], s, m, a], \[Beta]freq[n, \[Omega], Alm, s, m, a], {n, nInv+1}]
];

Leaver31Ang[\[Omega]_?NumericQ, Alm_?NumericQ, s_?IntegerQ, m_?IntegerQ, a_?NumericQ, nInv_?IntegerQ] := Module[{n},
  \[Beta]ang[nInv, \[Omega], Alm, s, m, a] + ContinuedFractionK[-\[Alpha]ang[nInv-n, \[Omega], s, m, a] \[Gamma]ang[nInv-n+1, \[Omega], s, m, a], \[Beta]ang[nInv-n,\[Omega], Alm, s, m, a], {n, 1, nInv}] + CF[-\[Alpha]ang[n-1, \[Omega], s, m, a] \[Gamma]ang[n, \[Omega], s, m, a], \[Beta]ang[n, \[Omega], Alm, s, m, a], {n, nInv+1}]
];


Options[QNMFrequencyLeaver] = {"InitialGuess" -> Automatic};


QNMFrequencyLeaver[s_Integer, l_Integer, m_Integer, n_, a_, OptionsPattern[]] :=
 Module[{\[Omega]guess, \[Omega]QNM, inInc, nInv = n, k, funcRad, funcAng, prec = prec[a]},
  \[Omega]guess=OptionValue["InitialGuess"];
  If[\[Omega]guess === Automatic,
    \[Omega]guess = Quiet[QNMFrequencyInterpolation[s, l, m, n][a], QNMFrequency::nointerp];
  ];
  If[\[Omega]guess == $Failed, 
    Message[QNMFrequency::nointerp, s, l, m, n];
    If[a==0,
      Which[
        l <= Max[n,2], 
        \[Omega]guess = SpectralInitialGuess[s, l, n];,
        l>n,
        \[Omega]guess = Schwarzfinit1[s, l, n];
      ];,
      \[Omega]guess = Kerrfinit[s, l, m, n, a];
    ];
  ];
  If[a==0,
    funcRad[\[Omega]_?NumericQ] := funcRad[\[Omega]] = Leaver[\[Omega], s, l, nInv];,
    funcRad[\[Omega]_?NumericQ] := funcRad[\[Omega]] = Leaver31[\[Omega], s, l, m, a, nInv];    
  ];
  \[Omega]QNM /. Check[FindRoot[funcRad[\[Omega]QNM]==0, {\[Omega]QNM, SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec], \[Omega]QNM -> $Failed, FindRoot::nlnum];
  \[Omega]QNM /. Quiet[Check[FindRoot[funcRad[\[Omega]QNM]==0, {\[Omega]QNM, SetPrecision[\[Omega]guess, prec]}, WorkingPrecision -> prec], \[Omega]QNM -> $Failed, FindRoot::nlnum], FindRoot::nlnum]
];


(* ::Subsection::Closed:: *)
(*QNMFrequency*)


SyntaxInformation[QNMFrequency] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[QNMFrequency] = {Method -> Automatic};


SetAttributes[QNMFrequency, {NumericFunction, Listable, NHoldAll}];


QNMFrequency[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] || !IntegerQ[n] || n < 0 :=
 (Message[QNMFrequency::params, s, l, m, n]; $Failed);


QNMFrequency[s_, l_, m_, n_, a_Complex, OptionsPattern[]] :=
 (Message[QNMFrequency::cmplx, a]; $Failed);


QNMFrequency[s_, l_, m_, n_, a_?InexactNumberQ, OptionsPattern[]] := 
 Module[{opts, \[Omega]ini, \[Omega]},
  Switch[OptionValue[Method],
    "Interpolation",
      \[Omega] = QNMFrequencyInterpolation[s, l, m, n][a],
    {"Interpolation", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyInterpolation]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyInterpolation[s, l, m, n, opts][a];,
    "IncidenceAmplitude",
      \[Omega] = QNMFrequencyInIncidenceAmplitude[s, l, m, n, a],
    {"IncidenceAmplitude", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyInIncidenceAmplitude]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyInIncidenceAmplitude[s, l, m, n, a, opts];,
    "Leaver",
      \[Omega] = QNMFrequencyLeaver[s, l, m, n, a],
    {"Leaver", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyLeaver]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyLeaver[s, l, m, n, a, opts];,
   Automatic | "SpheroidalEigenvalue",
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a];,
    {"SpheroidalEigenvalue", Rule[_,_]...},
	  opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyHyperboloidal]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyHyperboloidal[s, l, m, n, a, opts];,
    "Spectral1D",
      If[a!=0,
        Message[QNMFrequency::nokerr, "Spectral1D", a];
        \[Omega] = $Failed;
        ,
        \[Omega] = SpectralInitialGuess[s, l, n];
      ];,
    "Spectral2D",
      If[l =!= All || n =!= All,
        Message[QNMFrequency::allln, "Spectral2D", l, n];
        \[Omega] = $Failed;,
        \[Omega] = SpectralInitialGuessKerr[s, m, a];
      ];,
    {"Spectral2D", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[SpectralInitialGuessKerr]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      ];
      If[l!= All || n!= All,
        Message[QNMFrequency::allln, "Spectral2D", l, n];
        \[Omega] = $Failed;,
  	  \[Omega] = SpectralInitialGuessKerr[s, m, a, opts];
  	];,
    "Large-l Asymptotic",
      If[a == 0,
        \[Omega] = Schwarzfinit1[s, l, n];,
        \[Omega] = Kerrfinit[s, l, m, n, a];
      ];,
    "Large-n Asymptotic",
      If[a!=0,
        Message[QNMFrequency::nokerr, "Large-n Asymptotic", a];
        \[Omega] = $Failed;
        ,
        \[Omega] = Schwarzfinit2[n];
      ];,
    _,
      Message[QNMFrequency::optx, Method -> OptionValue[Method]];
      \[Omega] = $Failed;
  ];
  \[Omega]
];


QNMFrequency[s_, l_, m_, n_, 0, opts___] := QNMFrequency[s, l, m, n, 0.0, opts];


QNMFrequency /: N[QNMFrequency[s_, l_, m_, n_, a_?NumericQ, opts:OptionsPattern[]], Nopts___] :=
  QNMFrequency[s, l, m, n, N[a, Nopts], opts];


(* ::Section::Closed:: *)
(*Kerr-Newman: Hintz's separated system*)


(* ::Text:: *)
(*Quasinormal modes of Kerr-Newman black holes from P. Hintz's exact separation of the coupled gravito-electromagnetic perturbations (arXiv:2609.33661): an angular eigenvalue problem for the separation constant \[Lambda](a \[Omega]) and a first-order radial system for (w, w\[Sharp]).*)


(* ::Text:: *)
(*The quasinormal-mode boundary conditions in Hintz's variables, the radial solver (Frobenius series at the horizon, outgoing solution integrated inward along a complex ray, normalised Wronskian) and the spectral angular solver are adapted from A. Spiers, "Quasinormal modes from Hintz's separated Kerr-Newman equations: Reissner-Nordstr\[ODoubleDot]m tests and rotating black holes" (2026), its notebook kn_qnm_hintz.nb and spinning_KN/kn_angular.wl (functions hintzW, knSwY, knSpecMats, knSpecMatrix, knLamSpectral, knQNM and knTrack). Changes from that code: the horizon series is generated by a recurrence, to adaptive order and with the factor Exp[I \[Xi]0 r] removed (which makes it usable near Q = 1); the separation constant is split as \[Lambda] = b q + \[Delta] to avoid cancellations at small Q; the Wronskian normalisation is scale invariant; the ray angle is capped for strongly damped modes; an explicit Runge-Kutta method replaces NDSolve's default in machine precision (more accurate and faster); and the root search, the adaptive continuation from the Kerr mode and the precision handling are new.*)


(* ::Text:: *)
(*Units M = 1, modes ~ Exp[-I \[Omega] t + I m \[Phi]], a^2 + Q^2 < 1, Q > 0 (the spectrum is invariant under Q -> -Q). sg = +1 (-1) labels Hintz's positive (negative) spin system; the frequencies do not depend on it. b = +1 (-1) labels the gravitational-led (electromagnetic-led) family, which reduces to the s = \[PlusMinus]2 (s = \[PlusMinus]1) Teukolsky mode as Q -> 0.*)


(* ::Subsection::Closed:: *)
(*Angular problem*)


(* ::Text:: *)
(*On ker(S - i) the angular problem [Hintz, eq. (10)] reads I (\[PartialD]\[Theta] + \[CapitalTheta]0) u- + \[CapitalTheta]1 u+ = \[Lambda] u+, I (\[PartialD]\[Theta] - \[CapitalTheta]0) u+ - \[CapitalTheta]1 u- = \[Lambda] u-, with \[CapitalTheta]0 = 3/2 Cot[\[Theta]] + mm/Sin[\[Theta]] - c0 Sin[\[Theta]] and \[CapitalTheta]1 = -q - c1 Cos[\[Theta]]. The effective parameters are mm = sg m, c0 = sg a\[Omega] (1 + Q^2)/(1 - Q^2), c1 = sg 2 a\[Omega] Q/(1 - Q^2) and q = 3/(2Q).*)


KNEffectiveParameters[a\[Omega]_, Q_, m_, sg_] := {sg m, sg a\[Omega] (1 + Q^2)/(1 - Q^2), sg 2 a\[Omega] Q/(1 - Q^2), 3/(2 Q)};


(* ::Text:: *)
(*At a\[Omega] = 0 the eigenvalues are known exactly: \[Lambda] = b Sqrt[(l-1)(l+2) + q^2] for l >= 2 and \[Lambda] = -q for l = 1 (EM-led only). We work with \[Delta] = \[Lambda] - b q, which is small compared with q when Q is small.*)


KNDeltaRN[l_, Q_, b_] := If[l == 1, 0, b (l - 1) (l + 2)/(Sqrt[(l - 1) (l + 2) + 9/(4 Q^2)] + 3/(2 Q))];


(* ::Text:: *)
(*Spectral (Galerkin) method, from knSpecMats and knSpecMatrix of kn_angular.wl: u+/Sqrt[Sin[\[Theta]]] is expanded in normalised spin-weight-1 harmonics and u-/Sqrt[Sin[\[Theta]]] in spin-weight-2 harmonics with l <= L. The matrix is block tridiagonal in l. The matrix elements are computed by Gauss-Legendre quadrature with at least 60 digits, because evaluating high-degree Jacobi polynomials in machine precision cancels badly, and are then rounded to the working precision.*)


KNSwY[s_, l_, m_, x_] :=
 Module[{al = Abs[m + s], be = Abs[m - s], n},
  n = l - (al + be)/2;
  Sqrt[(2 n + al + be + 1) n! Gamma[n + al + be + 1]/(2^(al + be + 1) Gamma[n + al + 1] Gamma[n + be + 1])] *
    (1 - x)^(al/2) (1 + x)^(be/2) JacobiP[n, al, be, x]
];


KNSpectralMatrices[mm_, L_, wp_] := KNSpectralMatrices[mm, L, wp] =
 Module[{l1 = Max[1, Abs[mm]], l2 = Max[2, Abs[mm]], K = 2 L + 12, qp, xs, ws, x, f1, f2, g12, n1, n2, Mu12},
  qp = If[wp === MachinePrecision, 60, Max[60, Ceiling[2 wp]]];
  {xs, ws} = With[{gl = NIntegrate`GaussRuleData[K, qp]}, {2 gl[[1]] - 1, 2 gl[[2]]}];
  f1 = Table[Evaluate[Table[KNSwY[1, l, mm, x], {l, l1, L}] /. x -> xx], {xx, xs}];
  (* Sin[\[Theta]] times the spin-weight-2 basis functions *)
  g12 = Table[Evaluate[Table[((1 - x^2) D[KNSwY[1, l, mm, x], x] + (x + mm) KNSwY[1, l, mm, x])/Sqrt[(l - 1) (l + 2)],
      {l, l2, L}] /. x -> xx], {xx, xs}];
  f2 = g12/Sqrt[1 - xs^2];
  n1 = L - l1 + 1;
  n2 = L - l2 + 1;
  Mu12 = Table[If[l1 + i == l2 + j, Sqrt[(l2 + j - 2) (l2 + j + 1)], 0], {i, n1}, {j, n2}];
  N[{Mu12, Transpose[f1] . (ws xs f1), Transpose[f2] . (ws xs f2), Transpose[f1] . (ws g12)}, wp]
];


(* ::Text:: *)
(*The matrix is {{-q + A11, A12}, {A21, q + A22}}; its eigenvalues are \[Lambda].*)


KNSpectralBlocks[a\[Omega]_, Q_, m_, sg_, L_, wp_] :=
 Module[{mm, c0, c1, q, Mu12, C11, C22, S12},
  {mm, c0, c1, q} = KNEffectiveParameters[a\[Omega], Q, m, sg];
  {Mu12, C11, C22, S12} = KNSpectralMatrices[mm, L, wp];
  {-c1 C11, I Mu12 - I c0 S12, -I Transpose[Mu12] + I c0 Transpose[S12], c1 C22}
];


(* All values of \[Delta] = \[Lambda] - b q of the truncated problem *)
KNDeltaCandidates[{A11_, A12_, A21_, A22_}, q_, b_] :=
  Eigenvalues[ArrayFlatten[{{A11 - (1 + b) q IdentityMatrix[Length[A11]], A12}, {A21, A22 + (1 - b) q IdentityMatrix[Length[A22]]}}]];


(* ::Text:: *)
(*For |\[Delta]| << q the eigenvalue solver determines \[Delta] only to an absolute accuracy ~ q times the working precision. In machine precision we then refine it with the Schur complement, which does not contain q on the diagonal: \[Delta] is an eigenvalue of T(\[Delta]) = A22 + A21 . ((2q + \[Delta]) - A11)^-1 . A12 (b = +1) or T(\[Delta]) = A11 - A12 . ((2q - \[Delta]) + A22)^-1 . A21 (b = -1). The fixed-point iteration converges at a rate ~ |\[Delta]|/(2q).*)


KNDeltaRefine[{A11_, A12_, A21_, A22_}, q_, b_, \[Delta]0_, tol_] :=
 Module[{T, \[Delta] = \[Delta]0, \[Delta]new, change, prev = Infinity},
  T[d_] := If[b == 1,
    A22 + A21 . LinearSolve[(2 q + d) IdentityMatrix[Length[A11]] - A11, A12],
    A11 - A12 . LinearSolve[(2 q - d) IdentityMatrix[Length[A22]] + A22, A21]];
  Do[
    \[Delta]new = First[Nearest[Eigenvalues[T[\[Delta]]], \[Delta]]];
    change = Abs[\[Delta]new - \[Delta]];
    If[change > prev, Return[\[Delta]0, Module]];
    \[Delta] = \[Delta]new;
    If[change <= tol (Abs[\[Delta]] + 10^-3), Break[]];
    prev = change;
  , {8}];
  \[Delta]
];


(* \[Delta] on the branch closest to \[Delta]ref. The refinement is only needed (and only used) in machine precision: in
   extended precision the guard digits absorb the cancellation, and significance arithmetic in the Schur complement
   would underestimate the precision of the result. *)
KNDelta[a\[Omega]_, Q_, m_, sg_, b_, \[Delta]ref_, L_, wp_] :=
 Module[{blocks, q, \[Delta]},
  blocks = KNSpectralBlocks[a\[Omega], Q, m, sg, L, wp];
  q = 3/(2 Q);
  \[Delta] = First[Nearest[KNDeltaCandidates[blocks, q, b], \[Delta]ref]];
  If[wp === MachinePrecision,
    If[Abs[q] > 10 Abs[\[Delta]], \[Delta] = KNDeltaRefine[blocks, q, b, \[Delta], 10^(1 - $MachinePrecision)]],
    \[Delta] = SetPrecision[\[Delta], wp]
  ];
  \[Delta]
];


(* ::Text:: *)
(*Branch labels (l, m, b) are defined by continuity from a\[Omega] = 0, where the eigenvalue is known exactly, along the straight line to a\[Omega] (as lamBranch in kn_qnm_hintz.nb). This is done in machine precision; the result serves as the reference for picking the eigenvalue at nearby frequencies.*)


KNDeltaBranch[a\[Omega]_, Q_, m_, l_, b_, sg_, L_] :=
 Module[{\[Delta], \[Delta]prev, \[Delta]new, nst},
  \[Delta] = N[KNDeltaRN[l, Q, b]];
  If[a\[Omega] == 0, Return[\[Delta], Module]];
  nst = Max[10, Ceiling[Abs[a\[Omega]]/0.02]];
  \[Delta]prev = \[Delta];
  Do[
    \[Delta]new = KNDelta[k/nst a\[Omega], Q, m, sg, b, If[k == 1, \[Delta], 2 \[Delta] - \[Delta]prev], L, MachinePrecision];
    \[Delta]prev = \[Delta];
    \[Delta] = \[Delta]new;
  , {k, 1, nst}];
  \[Delta]
];


(* ::Subsection::Closed:: *)
(*Radial problem*)


(* ::Text:: *)
(*Hintz's radial system [eqs. (14)-(15) of arXiv:2609.33661]: \[PartialD]r w = P w\[Sharp], \[CapitalDelta] \[PartialD]r w\[Sharp] = Mm w - c w\[Sharp], with \[CapitalDelta] = r^2 - 2r + a^2 + Q^2, \[Xi]0 = 2\[Omega]Q^2/(1-Q^2), R1 = \[Rho]1(r) + I q, \[Rho]1 = -sg 2\[Omega]Q(r - Q^2)/(1-Q^2), P = R1 + I sg \[Lambda] = \[Rho]1 + I \[Lambda]p, Mm = R1 - I sg \[Lambda] = \[Rho]1 + I \[Lambda]m, \[Lambda]p = q + sg \[Lambda], \[Lambda]m = q - sg \[Lambda], c = (1 - 3sg)/2 \[CapitalDelta]' + 2I (K - \[Xi]0 \[CapitalDelta]) and K = -(r^2+a^2)\[Omega] + a m. The quasinormal-mode solution is smooth at r+ and behaves as Exp[I k r] r^p (1 + O(1/r)), k = 2\[Omega] + \[Xi]0, at infinity.*)


KNRadialParameters[\[Omega]_, a_, Q_, m_, sg_, b_, \[Delta]_] :=
 Module[{d, q = 3/(2 Q), \[Lambda]p, \[Lambda]m, \[Xi]0},
  d = 2 Sqrt[1 - a^2 - Q^2];
  (* \[Lambda] = b q + \[Delta], without cancellation in q \[PlusMinus] sg \[Lambda] *)
  \[Lambda]p = (1 + sg b) q + sg \[Delta];
  \[Lambda]m = (1 - sg b) q - sg \[Delta];
  \[Xi]0 = 2 \[Omega] Q^2/(1 - Q^2);
  <|"\[Omega]" -> \[Omega], "a" -> a, "Q" -> Q, "m" -> m, "sg" -> sg, "b" -> b, "\[Delta]" -> \[Delta], "\[Lambda]" -> b q + \[Delta],
    "\[Lambda]p" -> \[Lambda]p, "\[Lambda]m" -> \[Lambda]m, "d" -> d, "rp" -> 1 + d/2, "rm" -> 1 - d/2, "\[Xi]0" -> \[Xi]0,
    "\[Rho]1p" -> -sg 2 \[Omega] Q/(1 - Q^2), "k" -> 2 \[Omega] + \[Xi]0|>
];


(* Coefficients P, Mm, \[CapitalDelta], c of the radial system as functions of r *)
KNRadialCoefficients[p_Association] :=
 Module[{\[Omega], a, Q, m, sg, \[Lambda]p, \[Lambda]m, \[Xi]0, \[Rho]1p},
  {\[Omega], a, Q, m, sg, \[Lambda]p, \[Lambda]m, \[Xi]0, \[Rho]1p} = Lookup[p, {"\[Omega]", "a", "Q", "m", "sg", "\[Lambda]p", "\[Lambda]m", "\[Xi]0", "\[Rho]1p"}];
  {
   Function[r, \[Rho]1p (r - Q^2) + I \[Lambda]p],
   Function[r, \[Rho]1p (r - Q^2) + I \[Lambda]m],
   Function[r, r^2 - 2 r + a^2 + Q^2],
   Function[r, (1 - 3 sg) (r - 1) + 2 I (-(r^2 + a^2) \[Omega] + a m - \[Xi]0 (r^2 - 2 r + a^2 + Q^2))]
  }
];


(* ::Text:: *)
(*Integration method. The note used NDSolve's default method (LSODA) in machine precision; an explicit Runge-Kutta method is faster and, for the same precision goal, more accurate for this system (about 2*10^-12 -> 7*10^-14 at Q = 0.5 and 10^-8 -> 4*10^-10 at Q = 0.999 for the Wronskian zero). Its stiffness test is switched off: near Q = 1 it can stop the integration although the explicit method copes. In arbitrary precision the default (Adams) can fail at the first step for working precisions above about 40 digits; extrapolation is robust and fast there.*)


KNODEMethod[wp_] := If[wp === MachinePrecision, {"ExplicitRungeKutta", "StiffnessTest" -> False}, "Extrapolation"];


(* ::Text:: *)
(*Solution smooth at the horizon. Hintz's w contains the factor Exp[I \[Xi]0 r] [Hintz, eq. (14)], and \[Xi]0 = 2\[Omega]Q^2/(1-Q^2) is large near Q = 1; we therefore work with w~ = Exp[-I \[Xi]0 (r - r+)] w and w\[Sharp]~ = Exp[-I \[Xi]0 (r - r+)] w\[Sharp], which satisfy \[PartialD]r w~ = -I \[Xi]0 w~ + P w\[Sharp]~ and \[CapitalDelta] \[PartialD]r w\[Sharp]~ = Mm w~ - (c + I \[Xi]0 \[CapitalDelta]) w\[Sharp]~ and vary on the scale of the physical solution. We expand w~ = Sum[a_n x^n], w\[Sharp]~ = Sum[b_n x^n] with x = (r - r+)/(r+ - r-) and a_0 = 1. The only other finite singular point is r- (x = -1), so the series converges for |x| < 1. Substituting into the radial system gives the recurrence below (for \[Xi]0 = 0 it is equivalent to the linear system solved in hintzW). With nser = Automatic, terms are added until the series has converged at x = 1/3 to the working precision.*)


KNHorizonSeriesCoefficients[p_Association, nser_, wp_] :=
 Module[{d, rp, a, Q, m, sg, \[Omega], \[Xi]0, \[Rho]1p, \[Lambda]p, \[Lambda]m, \[Alpha]0, \[Alpha]1, \[Gamma]0, \[Gamma]1, \[Gamma]2,
         an, bn, n, nmax, eps, w, ws, xn, ta, tb, nsmall = 0, aprev, bprev, bprev2, anew, bnew, x1 = SetPrecision[1/3, wp]},
  {d, rp, a, Q, m, sg, \[Omega], \[Xi]0, \[Rho]1p, \[Lambda]p, \[Lambda]m} =
    Lookup[p, {"d", "rp", "a", "Q", "m", "sg", "\[Omega]", "\[Xi]0", "\[Rho]1p", "\[Lambda]p", "\[Lambda]m"}];
  (* P = \[Alpha]0 + I \[Lambda]p + \[Alpha]1 x, Mm = \[Alpha]0 + I \[Lambda]m + \[Alpha]1 x, c + I \[Xi]0 \[CapitalDelta] = \[Gamma]0 + \[Gamma]1 x + \[Gamma]2 x^2 *)
  \[Alpha]0 = \[Rho]1p (rp - Q^2);
  \[Alpha]1 = d \[Rho]1p;
  \[Gamma]0 = (1 - 3 sg)/2 d + 2 I (a m - (rp^2 + a^2) \[Omega]);
  \[Gamma]1 = (1 - 3 sg) d - I d (4 rp \[Omega] + \[Xi]0 d);
  \[Gamma]2 = -I d^2 (2 \[Omega] + \[Xi]0);
  nmax = If[IntegerQ[nser], nser, 2000];
  eps = 10^-If[wp === MachinePrecision, $MachinePrecision, wp];
  {aprev, bprev, bprev2} = {1, (\[Alpha]0 + I \[Lambda]m)/\[Gamma]0, 0};
  an = Internal`Bag[{aprev}];
  bn = Internal`Bag[{bprev}];
  {w, ws} = {aprev, bprev};
  xn = 1;
  For[n = 1, n <= nmax, n++,
    anew = (d ((\[Alpha]0 + I \[Lambda]p) bprev + \[Alpha]1 bprev2) - I \[Xi]0 d aprev)/n;
    bnew = ((\[Alpha]0 + I \[Lambda]m) anew + \[Alpha]1 aprev - (d (n - 1) + \[Gamma]1) bprev - \[Gamma]2 bprev2)/(d n + \[Gamma]0);
    Internal`StuffBag[an, anew];
    Internal`StuffBag[bn, bnew];
    {aprev, bprev2, bprev} = {anew, bprev, bnew};
    If[!IntegerQ[nser],
      xn *= x1;
      ta = anew xn;
      tb = bnew xn;
      w += ta;
      ws += tb;
      If[Abs[ta] <= eps Abs[w] && Abs[tb] <= eps Abs[ws], nsmall++, nsmall = 0];
      If[nsmall >= 4 && n >= 8, Break[]];
    ];
  ];
  {Internal`BagPart[an, All], Internal`BagPart[bn, All]}
];


(* The series and its largest term relative to its sum (a measure of cancellation) at x0 *)
KNHorizonSeriesValue[{A_, B_}, x0_] :=
 Module[{xs = x0^Range[0, Length[A] - 1], w, ws},
  w = A . xs;
  ws = B . xs;
  {w, ws, Max[Max[Abs[A xs]]/Abs[w], If[ws == 0, 0, Max[Abs[B xs]]/Abs[ws]]]}
];


(* ::Text:: *)
(*Horizon solution: the series is evaluated at x0 and the system for (w, w\[Sharp]) is integrated along the real axis to the matching radius rmatch. By default x0 = 1/3, as in the note. If the series has large intermediate terms at x = 1/3 (its sum would cancel), x0 is reduced by factors of 3 until the largest term is below 10^6 times the sum (10^(wp - pg) in extended precision); with the factor Exp[I \[Xi]0 r] removed this is rarely needed. Starting closer to the horizon is not better in general, because the integration from there loses accuracy. The result contains the series coefficients of w~ and w\[Sharp]~ (valid for r+ <= r <= r0, with w = Exp[I \[Xi]0 (r - r+)] w~) and the solution (w, w\[Sharp]) on [r0, rmatch].*)


KNHorizonSolution[p_Association, nser_, x00_, rmatch_, wp_, pg_] :=
 Module[{coeffs, x0, ser, cancel, P, Mm, \[CapitalDelta], c, ef, r0, w, ws, r, sol},
  coeffs = KNHorizonSeriesCoefficients[p, nser, wp];
  If[x00 === Automatic,
    cancel = 10^Max[6, If[wp === MachinePrecision, 0, wp - pg]];
    x0 = SetPrecision[1/3, wp];
    ser = KNHorizonSeriesValue[coeffs, x0];
    While[ser[[3]] > cancel && x0 > 1/1000,
      x0 = x0/3;
      ser = KNHorizonSeriesValue[coeffs, x0];
    ];
    ,
    x0 = SetPrecision[x00, wp];
    ser = KNHorizonSeriesValue[coeffs, x0];
  ];
  {P, Mm, \[CapitalDelta], c} = KNRadialCoefficients[p];
  r0 = p["rp"] + p["d"] x0;
  (* w = Exp[I \[Xi]0 (r - r+)] w~; the system for w itself is integrated (the one for w~ is stiff near extremality) *)
  ef = Exp[I p["\[Xi]0"] p["d"] x0];
  sol = Quiet[First @ NDSolve[{
      w'[r] == P[r] ws[r],
      ws'[r] == (Mm[r] w[r] - c[r] ws[r])/\[CapitalDelta][r],
      w[r0] == SetPrecision[ef ser[[1]], wp], ws[r0] == SetPrecision[ef ser[[2]], wp]},
    {w, ws}, {r, r0, rmatch},
    WorkingPrecision -> wp, PrecisionGoal -> pg, AccuracyGoal -> Infinity, MaxSteps -> Infinity, Method -> KNODEMethod[wp]],
    {NDSolve::ndsz, NDSolve::precw}];
  <|"SeriesCoefficients" -> coeffs, "x0" -> x0, "r0" -> r0, "rmatch" -> rmatch, "w" -> (w /. sol), "w\[Sharp]" -> (ws /. sol)|>
];


(* ::Text:: *)
(*Outgoing solution: integrated inward along the ray r = rmatch + s Exp[I \[Phi]] from s = \[Rho] to s = 0, where the outgoing solution is dominant. We take \[Phi] = \[Pi]/2 - Arg[\[Omega]], capped at 3\[Pi]/4 (5\[Pi]/4 for Re[\[Omega]] < 0) so that the ray stays away from r\[PlusMinus] for strongly damped modes. On the ray the exponential is factored out, w = Exp[I k r] u, w\[Sharp] = Exp[I k r] u\[Sharp], and the ray is parametrised by v = \[Rho] - s. The starting values are the leading asymptotic behaviour; the admixture of the ingoing solution is suppressed by Exp[-2 |\[Omega]| \[Rho] Cos[\[Phi] - \[Phi]opt]].*)


KNRayAngle[\[Omega]_] :=
 Module[{\[Phi] = \[Pi]/2 - Arg[\[Omega]]},
  If[Re[\[Omega]] >= 0, Min[\[Phi], 3 \[Pi]/4], Max[\[Phi], 5 \[Pi]/4]]
];


KNRayLength[\[Omega]_, digits_] :=
  Max[50, (digits Log[10] + 3)/(2 Abs[\[Omega]] Cos[KNRayAngle[\[Omega]] - (\[Pi]/2 - Arg[\[Omega]])])];


KNInfinitySolution[p_Association, rmatch_, \[Rho]_, wp_, pg_] :=
 Module[{P, Mm, \[CapitalDelta], c, \[Phi], e, k, rinf, ray, u, us, v, sol},
  {P, Mm, \[CapitalDelta], c} = KNRadialCoefficients[p];
  k = p["k"];
  \[Phi] = SetPrecision[KNRayAngle[p["\[Omega]"]], wp];
  e = Exp[I \[Phi]];
  rinf = rmatch + \[Rho] e;
  ray[vv_] := rmatch + (\[Rho] - vv) e;
  sol = Quiet[First @ NDSolve[{
      u'[v] == -e (P[ray[v]] us[v] - I k u[v]),
      us'[v] == -e (Mm[ray[v]] u[v] - (c[ray[v]] + I k \[CapitalDelta][ray[v]]) us[v])/\[CapitalDelta][ray[v]],
      u[0] == 1, us[0] == I k/P[rinf]},
    {u, us}, {v, 0, \[Rho]},
    WorkingPrecision -> wp, PrecisionGoal -> pg, AccuracyGoal -> Infinity, MaxSteps -> Infinity, Method -> KNODEMethod[wp]],
    {NDSolve::ndsz, NDSolve::precw}];
  <|"rmatch" -> rmatch, "\[Rho]" -> \[Rho], "\[Phi]" -> \[Phi], "k" -> k, "u" -> (u /. sol), "u\[Sharp]" -> (us /. sol)|>
];


(* ::Text:: *)
(*Both solutions, packaged for reuse (e.g. for constructing radial functions): the horizon solution (series and real-axis solution up to rmatch) and the outgoing solution on the ray, (u, u\[Sharp]) = Exp[-I k r] (w, w\[Sharp]) as functions of v, with r = rmatch + (\[Rho] - v) Exp[I \[Phi]]. At a quasinormal frequency the two are proportional at rmatch.*)


Options[KNRadialSolutions] = {"SeriesOrder" -> Automatic, "SeriesPoint" -> Automatic, "MatchingRadius" -> 6, "RayLength" -> Automatic,
  WorkingPrecision -> MachinePrecision, PrecisionGoal -> Automatic};


KNRadialSolutions[p_Association, OptionsPattern[]] :=
 Module[{wp, pg, rmatch, \[Rho], hor, inf},
  wp = OptionValue[WorkingPrecision];
  pg = OptionValue[PrecisionGoal];
  If[pg === Automatic, pg = If[wp === MachinePrecision, 11, wp - 8]];
  rmatch = SetPrecision[OptionValue["MatchingRadius"], wp];
  \[Rho] = OptionValue["RayLength"];
  If[\[Rho] === Automatic, \[Rho] = KNRayLength[p["\[Omega]"], pg + 1]];
  \[Rho] = SetPrecision[\[Rho], wp];
  hor = KNHorizonSolution[p, OptionValue["SeriesOrder"], OptionValue["SeriesPoint"], rmatch, wp, pg];
  inf = KNInfinitySolution[p, rmatch, \[Rho], wp, pg];
  <|"Parameters" -> p, "Horizon" -> hor, "Infinity" -> inf|>
];


(* ::Text:: *)
(*Normalised Wronskian of the two solutions at rmatch. Its zeros in \[Omega] are the quasinormal modes; the normalisation makes it invariant under rescaling of w\[Sharp].*)


KNWronskian[sols_Association] :=
 Module[{rmatch, \[Rho], wH, wsH, uI, usI},
  rmatch = sols["Horizon"]["rmatch"];
  \[Rho] = sols["Infinity"]["\[Rho]"];
  (* an integration that stopped early (e.g. for a wild iterate of the root search) gives no value *)
  If[!TrueQ[Abs[sols["Horizon"]["w"]["Domain"][[1, 2]] - rmatch] <= 10^-6 rmatch] ||
     !TrueQ[Abs[sols["Infinity"]["u"]["Domain"][[1, 2]] - \[Rho]] <= 10^-6 \[Rho]],
    Return[$Failed, Module]];
  (* the factor Exp[I k rmatch] of the ray solution is common to w and w\[Sharp] and drops out *)
  wH = sols["Horizon"]["w"][rmatch];
  wsH = sols["Horizon"]["w\[Sharp]"][rmatch];
  uI = sols["Infinity"]["u"][\[Rho]];
  usI = sols["Infinity"]["u\[Sharp]"][\[Rho]];
  (wH usI - wsH uI)/(Abs[wH usI] + Abs[wsH uI])
];


(* ::Subsection::Closed:: *)
(*Root finding and continuation*)


(* ::Text:: *)
(*Secant iteration. The Wronskian has a noise floor set by the working precision and the precision goal, so besides the step criterion (step <= tol |z|) we stop after four consecutive small steps (< Sqrt[tol] |z|) that do not reach tol, and return the iterate with the smallest residual. Steps are limited to a tenth of |z|, and the iteration is abandoned when it moves further than maxdev from the starting value. Iterates are kept at the working precision wp. Returns {root, converged, last step}; the last step estimates the accuracy.*)


KNSecant[f_, z0_, z1_, tol_, maxit_, maxdev_, wp_] :=
 Module[{za = z0, zb = z1, fa, fb, zc, fc, best, fbest, small = 0, conv = False, dz, last = Infinity},
  fa = f[za];
  fb = f[zb];
  If[!NumericQ[fa] || !NumericQ[fb], Return[{$Failed, False, Infinity}, Module]];
  {best, fbest} = If[Abs[fa] <= Abs[fb], {za, fa}, {zb, fb}];
  Do[
    If[fb == fa, Break[]];
    dz = -fb (zb - za)/(fb - fa);
    If[Abs[dz] > Abs[zb]/10, dz = dz Abs[zb]/(10 Abs[dz])];
    zc = zb + dz;
    If[wp =!= MachinePrecision, zc = SetPrecision[zc, wp]];
    If[Abs[zc - z0] > maxdev, Break[]];
    fc = f[zc];
    If[$KNDebug, Print["    secant: ", zc, "  |W| = ", N[Abs[fc]]]];
    If[!NumericQ[fc], Break[]];
    last = Abs[zc - zb];
    If[Abs[fc] < Abs[fbest], {best, fbest} = {zc, fc}];
    If[last <= tol Abs[zc], {best, conv} = {zc, True}; Break[]];
    (* noise floor: the iterates stay within Sqrt[tol] |z| of each other but do not converge further *)
    If[last < Sqrt[tol] Abs[zc], small++, small = 0];
    If[small >= 4, conv = True; Break[]];
    {za, fa, zb, fb} = {zb, fb, zc, fc};
  , {maxit}];
  {best, conv, last}
];


(* ::Text:: *)
(*Root of the Wronskian at fixed (a, Q), starting from \[Omega]guess. The angular branch is fixed by continuation from a\[Omega] = 0 at \[Omega]guess, and the separation constant at nearby \[Omega] is the eigenvalue closest to it. Returns {\[Omega], \[Delta], converged, last step}.*)


Options[KNRoot] = {"SeriesOrder" -> Automatic, "MatchingRadius" -> 6, "RayLength" -> Automatic,
  "AngularTruncation" -> Automatic, WorkingPrecision -> MachinePrecision, PrecisionGoal -> Automatic,
  "Tolerance" -> Automatic, "MaxIterations" -> 30, "MaxDeviation" -> Infinity, "BranchReference" -> Automatic};


KNRoot[l_, m_, b_, sg_, a_, Q_, \[Omega]guess_, opts:OptionsPattern[]] :=
 Module[{wp, L, \[Delta]ref, \[Delta]of, W, radopts, tol, res, \[Omega]0, \[Omega]1},
  wp = OptionValue[WorkingPrecision];
  L = OptionValue["AngularTruncation"];
  If[L === Automatic, L = Max[l, Abs[m], 2] + 22];
  \[Delta]ref = OptionValue["BranchReference"];
  If[\[Delta]ref === Automatic, \[Delta]ref = KNDeltaBranch[SetPrecision[a \[Omega]guess, MachinePrecision], SetPrecision[Q, MachinePrecision], m, l, b, sg, L]];
  \[Delta]of[\[Omega]_] := \[Delta]of[\[Omega]] =
    If[a == 0, N[KNDeltaRN[l, Q, b], wp], KNDelta[a \[Omega], Q, m, sg, b, \[Delta]ref, L, wp]];
  radopts = FilterRules[{opts}, Options[KNRadialSolutions]];
  W[\[Omega]_?NumericQ] := W[\[Omega]] =
    KNWronskian[KNRadialSolutions[KNRadialParameters[\[Omega], a, Q, m, sg, b, \[Delta]of[\[Omega]]], radopts]];
  tol = OptionValue["Tolerance"];
  If[tol === Automatic, tol = If[wp === MachinePrecision, 10^-12, 10^(3 - wp)]];
  \[Omega]0 = SetPrecision[\[Omega]guess, wp];
  (* second starting point: a small fraction of the damping rate, which sets the scale on which W varies *)
  \[Omega]1 = \[Omega]0 + SetPrecision[If[wp === MachinePrecision, 10^-3, 10^-8] Abs[Im[\[Omega]0]], wp];
  res = KNSecant[W, \[Omega]0, \[Omega]1, tol, OptionValue["MaxIterations"], OptionValue["MaxDeviation"], wp];
  If[res[[1]] === $Failed, Return[{$Failed, $Failed, False, Infinity}, Module]];
  {res[[1]], \[Delta]of[res[[1]]], res[[2]], res[[3]]}
];


(* ::Text:: *)
(*Continuation in the charge at fixed spin, starting from the Kerr frequency at Q = 0 (computed with QNMFrequency). This follows the continuation used for the KN tables of the note (knTrack), with adaptive steps. The continuation variable is \[Tau] = Sqrt[1 - a^2 - Q^2] (half the horizon separation): the frequency is analytic in Q^2 near Q = 0 and, near extremality, varies smoothly with \[Tau] but not with Q. The predictor is the interpolating polynomial in \[Tau] through the last (at most three) computed points, including the Kerr point. A step is rejected and halved when the corrector moves the frequency by more than StepTolerance |Im[\[Omega]]| from the prediction (or does not converge), and the next step is doubled when the prediction was accurate; a step never reduces \[Tau] by more than half.*)


$KNDebug = False;


Options[KNContinuation] = Join[{"KerrFrequency" -> Automatic, "InitialCharge" -> 1/20, "StepTolerance" -> 1/10, "MinStep" -> 10^-4}, Options[KNRoot]];


KNContinuation[s_, l_, m_, n_, a_, Qt_, sg_, opts:OptionsPattern[]] :=
 Module[{b, \[Omega]K, Qs, pts, rootopts, solve, \[Tau], \[Tau]of, \[Tau]T, \[Tau]new, Qnew, h, pred, res, dev, \[Eta], hmin, last},
  b = If[Abs[s] == 2, 1, -1];
  rootopts = FilterRules[{opts}, Options[KNRoot]];
  \[Eta] = OptionValue["StepTolerance"];
  (* Kerr frequency; the sign of s does not affect the frequency and the s < 0 data are tabulated. It only serves as
     a starting value (the first root at small Q is checked against it), so messages from QNMFrequency are suppressed.
     Without tabulated data (l > 7 or n > 6) QNMFrequency starts from asymptotic approximations, which can select a
     different mode; the Kerr frequency can then be supplied with the KerrFrequency option. *)
  \[Omega]K = OptionValue["KerrFrequency"];
  If[\[Omega]K === Automatic, \[Omega]K = Quiet[QNMFrequency[-Abs[s], l, m, n, a]], \[Omega]K = SetPrecision[\[Omega]K, MachinePrecision]];
  If[!NumericQ[\[Omega]K], Message[QNMFrequencyKN::seed, -Abs[s], l, m, n, a]; Return[$Failed, Module]];
  solve[Q_, guess_] := KNRoot[l, m, b, sg, a, Q, guess, "MaxDeviation" -> 2 \[Eta] Abs[Im[guess]], "MaxIterations" -> 12, rootopts];
  \[Tau]of[Q_] := Sqrt[1 - a^2 - Q^2];
  (* first point at small charge, from the Kerr frequency (which differs from it by O(Q^2)). For high frequencies the
     root search converges only from close by, so Qs is halved until the root is found. *)
  Qs = Min[Qt, N[OptionValue["InitialCharge"]], 0.3 Sqrt[1 - a^2]];
  While[True,
    res = solve[Qs, \[Omega]K];
    If[$KNDebug, Print["  Q = ", Qs, ": ", res]];
    If[TrueQ[res[[3]]] && Abs[res[[1]] - \[Omega]K] <= \[Eta] Abs[Im[\[Omega]K]], Break[]];
    Qs = Qs/2;
    If[Qs < 10^-3, Message[QNMFrequencyKN::cont, a, 2 Qs]; Return[$Failed, Module]];
  ];
  pts = {{\[Tau]of[0], \[Omega]K}, {\[Tau]of[Qs], res[[1]]}};
  last = res[[4]];
  \[Tau] = \[Tau]of[Qs];
  \[Tau]T = \[Tau]of[Qt];
  hmin = OptionValue["MinStep"] \[Tau]T;
  h = \[Tau] - \[Tau]T;
  While[\[Tau] > \[Tau]T,
    (* at most halve \[Tau] in one step: near extremality the frequency varies on the scale \[Tau] *)
    h = Min[h, \[Tau]/2];
    \[Tau]new = Max[\[Tau] - h, \[Tau]T];
    Qnew = If[\[Tau]new == \[Tau]T, Qt, Sqrt[1 - a^2 - \[Tau]new^2]];
    pred = InterpolatingPolynomial[Take[pts, -Min[3, Length[pts]]], \[Tau]new];
    res = solve[Qnew, pred];
    dev = If[TrueQ[res[[3]]], Abs[res[[1]] - pred], Infinity];
    If[$KNDebug, Print["  Q = ", Qnew, ", \[Tau] = ", \[Tau]new, ": pred ", pred, ", root ", res[[1]], ", conv ", res[[3]], ", dev/|Im| = ", dev/Abs[Im[pred]]]];
    If[dev <= \[Eta] Abs[Im[res[[1]]]],
      AppendTo[pts, {\[Tau]new, res[[1]]}];
      last = res[[4]];
      \[Tau] = \[Tau]new;
      If[dev <= \[Eta]/8 Abs[Im[res[[1]]]], h = 2 h];
      ,
      h = (\[Tau] - \[Tau]new)/2;
      If[h < hmin, Message[QNMFrequencyKN::cont, a, Qnew]; Return[$Failed, Module]];
    ];
  ];
  {pts[[-1, 2]], last}
];


(* ::Subsection::Closed:: *)
(*Hintz separated method*)


Options[QNMFrequencyKNHintz] = {"InitialGuess" -> Automatic, "KerrFrequency" -> Automatic, "SeriesOrder" -> Automatic, "MatchingRadius" -> 6,
  "RayLength" -> Automatic, "AngularTruncation" -> Automatic, "SpinSystem" -> 1, "AccuracyCheck" -> Automatic,
  WorkingPrecision -> Automatic, PrecisionGoal -> Automatic};


KNPrecision[x_] := Which[!InexactNumberQ[x], Infinity, x == 0 && Precision[x] == MachinePrecision, MachinePrecision, x == 0, Accuracy[x], True, Precision[x]];


QNMFrequencyKNHintz[s_, l_, m0_, n_, a0_, Q0_, opts:OptionsPattern[]] :=
 Module[{prec, wp, pg, m, a, Q, aM, QM, b, sg, L, guess, rootopts, \[Omega], err, res, wpi},
  (* precision of the result follows the input; the internal precision can be raised with WorkingPrecision *)
  prec = Min[KNPrecision[a0], KNPrecision[Q0]];
  wp = OptionValue[WorkingPrecision];
  If[wp === Automatic, wp = prec];
  (* symmetries: (a, m) -> (-a, -m) and Q -> -Q leave the spectrum unchanged *)
  m = If[a0 < 0, -m0, m0];
  a = Abs[a0];
  Q = Abs[Q0];
  b = If[Abs[s] == 2, 1, -1];
  sg = OptionValue["SpinSystem"];
  L = OptionValue["AngularTruncation"];
  If[L === Automatic, L = Max[l, Abs[m], 2] + 22];
  rootopts = {"SeriesOrder" -> OptionValue["SeriesOrder"], "MatchingRadius" -> OptionValue["MatchingRadius"],
    "RayLength" -> OptionValue["RayLength"], "AngularTruncation" -> L};
  pg = OptionValue[PrecisionGoal];
  (* machine-precision stage: continuation from the Kerr mode, or a root search from a user-supplied guess *)
  {aM, QM} = SetPrecision[{a, Q}, MachinePrecision];
  guess = OptionValue["InitialGuess"];
  If[guess === Automatic,
    res = KNContinuation[s, l, m, n, aM, QM, sg, rootopts, "KerrFrequency" -> OptionValue["KerrFrequency"],
      PrecisionGoal -> If[pg === Automatic, Automatic, Min[pg, 13]]];
    If[res === $Failed, Return[$Failed, Module]];
    {\[Omega], err} = res;
    ,
    res = KNRoot[l, m, b, sg, aM, QM, SetPrecision[guess, MachinePrecision], rootopts, PrecisionGoal -> If[pg === Automatic, Automatic, Min[pg, 13]]];
    If[res[[1]] === $Failed || !TrueQ[res[[3]]], Message[QNMFrequencyKN::findroot, N[guess]]; Return[$Failed, Module]];
    {\[Omega], err} = res[[{1, 4}]];
  ];
  (* extended-precision stage: polish the root. NDSolve needs a working precision well above its precision goal. *)
  If[wp =!= MachinePrecision && wp > $MachinePrecision,
    wpi = Ceiling[wp] + 12;
    res = KNRoot[l, m, b, sg, SetPrecision[a, wpi], SetPrecision[Q, wpi], SetPrecision[\[Omega], wpi],
      Sequence @@ rootopts, WorkingPrecision -> wpi, PrecisionGoal -> If[pg === Automatic, Ceiling[wp] + 2, pg], "Tolerance" -> 10^-(Ceiling[wp] + 1),
      "BranchReference" -> KNDeltaBranch[aM \[Omega], QM, m, l, b, sg, L]];
    If[res[[1]] === $Failed, Message[QNMFrequencyKN::findroot, \[Omega]]; Return[$Failed, Module]];
    {\[Omega], err} = res[[{1, 4}]];
    If[!TrueQ[res[[3]]] || err > 10^-(wp - 2) Abs[\[Omega]], Message[QNMFrequencyKN::acc, N[\[Omega]], N[err], wp]];
    ,
    (* machine precision: near extremality the result can be much less accurate than ~10^-11, because the radial system
       is badly conditioned there, and the size of the last secant step does not show this. The error is then dominated
       by the precision goal of NDSolve, so we repeat the root search with the precision goal raised by one digit, keep
       that more accurate root, and use the change as an (upper) estimate of the error. *)
    If[TrueQ[OptionValue["AccuracyCheck"]] || (OptionValue["AccuracyCheck"] === Automatic && (1 - a^2 - Q^2 < 1/50 || Q > 19/20)),
      res = KNRoot[l, m, b, sg, aM, QM, \[Omega], rootopts, PrecisionGoal -> If[pg === Automatic, 12, Min[pg + 1, 13]]];
      If[NumericQ[res[[1]]] && TrueQ[res[[3]]],
        err = Abs[res[[1]] - \[Omega]];
        \[Omega] = res[[1]];
        ,
        err = Infinity;
      ];
      If[err > 10^-9 Abs[\[Omega]], Message[QNMFrequencyKN::acc, \[Omega], N[err], MachinePrecision]];
    ];
  ];
  If[prec === MachinePrecision, N[\[Omega]], SetPrecision[\[Omega], prec]]
];


(* ::Subsection::Closed:: *)
(*QNMFrequencyKN*)


SyntaxInformation[QNMFrequencyKN] =
 {"ArgumentsPattern" -> {_, _, _, _, _, _, OptionsPattern[]}};


Options[QNMFrequencyKN] = {Method -> Automatic};


SetAttributes[QNMFrequencyKN, {NumericFunction, Listable, NHoldAll}];


QNMFrequencyKN[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_, Q_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] || !IntegerQ[n] || n < 0 :=
 (Message[QNMFrequencyKN::params, s, l, m, n]; $Failed);


QNMFrequencyKN[s_?NumericQ, l_, m_, n_, a_, Q_, OptionsPattern[]] /; !MemberQ[{1, 2}, Abs[s]] :=
 (Message[QNMFrequencyKN::spin, s]; $Failed);


QNMFrequencyKN[s_, l_, m_, n_, a_, Q_, OptionsPattern[]] /; MatchQ[a, _Complex] || MatchQ[Q, _Complex] :=
 (Message[QNMFrequencyKN::cmplx, a, Q]; $Failed);


QNMFrequencyKN[s_, l_, m_, n_, a_?NumericQ, Q_?NumericQ, OptionsPattern[]] /;
  (InexactNumberQ[a] || InexactNumberQ[Q]) && TrueQ[a^2 + Q^2 >= 1] :=
 (Message[QNMFrequencyKN::extremal, a, Q]; $Failed);


QNMFrequencyKN[s_Integer, l_Integer, m_Integer, n_Integer, a_?NumericQ, Q_?NumericQ, OptionsPattern[]] /;
  (InexactNumberQ[a] || InexactNumberQ[Q]) && Q == 0 :=
 QNMFrequency[-Abs[s], l, m, n, If[InexactNumberQ[a], a, N[a, Min[KNPrecision[a], KNPrecision[Q]]]]];


QNMFrequencyKN[s_Integer, l_Integer, m_Integer, n_Integer, a_?NumericQ, Q_?NumericQ, OptionsPattern[]] /;
  (InexactNumberQ[a] || InexactNumberQ[Q]) :=
 Module[{opts, \[Omega]},
  Switch[OptionValue[Method],
    Automatic | "HintzSeparated",
      \[Omega] = QNMFrequencyKNHintz[s, l, m, n, a, Q];,
    {"HintzSeparated", Rule[_, _]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMFrequencyKNHintz]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMFrequencyKN::optx, Method -> OptionValue[Method]];
      ];
      \[Omega] = QNMFrequencyKNHintz[s, l, m, n, a, Q, opts];,
    _,
      Message[QNMFrequencyKN::optx, Method -> OptionValue[Method]];
      \[Omega] = $Failed;
  ];
  \[Omega]
];


QNMFrequencyKN[s_, l_, m_, n_, 0, 0, opts___] := QNMFrequencyKN[s, l, m, n, 0.0, 0.0, opts];


QNMFrequencyKN /: N[QNMFrequencyKN[s_, l_, m_, n_, a_?NumericQ, Q_?NumericQ, opts:OptionsPattern[]], Nopts___] :=
  QNMFrequencyKN[s, l, m, n, N[a, Nopts], N[Q, Nopts], opts];


(* ::Section::Closed:: *)
(*QNMRadial*)


(* ::Subsection::Closed:: *)
(*QNMRadialHyperboloidal*)


Options[QNMRadialHyperboloidal] = {
  "Coordinates" -> "Hyperboloidal",
  "NumPoints" -> 100
};


QNMRadialHyperboloidal[s_, l_, m_, n_, a_, \[Omega]_, opts:OptionsPattern[]] :=
 Module[{\[Lambda], ef, ns, Mat, RadialFunction, h, h\[Phi], numpoints, coords, domain},
  (* Load options values *)
  numpoints = OptionValue["NumPoints"];

  (* Check a valid Coordinates option has been specified *)
  coords = OptionValue["Coordinates"];
  If[!MemberQ[{"BL", "Boyer-Lindquist", "BoyerLindquist", "Hyperboloidal", "CompactifiedHyperboloidal"}, coords],
    Message[QNMRadial::coords, coords];
    Return[$Failed];
  ];

  (* Construct the discretized system *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  Mat = \[ScriptCapitalM][s, m, a, \[Omega], \[Lambda], numpoints];

  (* Calculate the eigenvectors *)
  ns = NullSpace[Mat];
  If[Length[ns]==0,
    ef = Eigensystem[Mat][[2, -1]];,
    ef = First[ns];
  ];

  Switch[coords,
  "Hyperboloidal",
    RadialFunction = Function[{r}, Evaluate[chebInterp[Reverse[ef/ef[[-1]]], {0, 1/rp[a, M]}][1/r]]];
    domain = {rp[a, M], \[Infinity]};,
  "CompactifiedHyperboloidal",
    RadialFunction = Function[{\[Sigma]}, Evaluate[chebInterp[Reverse[ef/ef[[-1]]], {0, 1/rp[a, M]}][\[Sigma]]]];
    domain = {0, 1/rp[a,M]};,
  "BL" | "BoyerLindquist" | "Boyer-Lindquist",
    RadialFunction = Function[{r}, Evaluate[
      With[{rp = rp[a, M], rm = rm[a, M]},
        h = (2 M rp )/(rp-rm) Log[r-rp]-(2 M rm )/(rp-rm) Log[r-rm]-r-4 M Log[r];
        h\[Phi] = a/(rp-rm) Log[(r-rp)/(r-rm)];
        chebInterp[Reverse[ef/ef[[-1]]], {0, 1/rp[a, M]}][1/r] Exp[-I*\[Omega]*h+I*m*h\[Phi]]]]];
    domain = {rp[a, M], \[Infinity]};,
  _,
    Message[QNMRadial::coords, coords];
    Return[$Failed];
  ];

  (* Return QNMRadialFunction *)
  QNMRadialFunction[<|"s" -> s, "l" -> l, "m" -> m, "n" -> n, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
    "Method" -> "SpectralHyperboloidal", "Amplitudes" -> <|"\[ScriptCapitalH]" -> ef[[-1]]/ef[[-1]], "\[ScriptCapitalI]" -> ef[[1]]/ef[[-1]]|>, "RadialFunction" -> RadialFunction,
    "Coordinates" -> coords, "Domain" -> domain|>]
]


(* ::Subsection::Closed:: *)
(*QNMRadial*)


SyntaxInformation[QNMRadial] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[QNMRadial] = {
  Method -> Automatic,
  "Frequency" -> Automatic
};


SetAttributes[QNMRadial, {Listable, NHoldAll}];


QNMRadial[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] || !IntegerQ[n] || n < 0 :=
 (Message[QNMRadial::params, s, l, m, n]; $Failed);


QNMRadial[s_, l_, m_, n_, a_Complex, OptionsPattern[]] :=
 (Message[QNMRadial::cmplx, a]; $Failed);


QNMRadial[s_?NumericQ, l_?NumericQ, m_?NumericQ, n_?NumericQ, a_?InexactNumberQ, OptionsPattern[]] :=
 Module[{\[Omega], opts, qnm},
  (* Get the frequency *)
  If[OptionValue["Frequency"] === Automatic,
    \[Omega] = QNMFrequency[s, l, m, n, a];,
    \[Omega] = OptionValue["Frequency"];
  ];

  Switch[OptionValue[Method],
    Automatic|"Hyperboloidal",
      qnm = QNMRadialHyperboloidal[s, l, m, n, a, \[Omega]],
    {"Hyperboloidal", Rule[_,_]...},
      opts = FilterRules[Rest[OptionValue[Method]], Options[QNMRadialHyperboloidal]];
      If[opts =!= Rest[OptionValue[Method]],
        Message[QNMRadial::optx, Method -> OptionValue[Method]];
      ];
      qnm = QNMRadialHyperboloidal[s, l, m, n, a, \[Omega], opts],
    _,
     Message[QNMRadial::optx, Method -> OptionValue[Method]];
     qnm = $Failed;
  ];
  qnm
]


QNMRadial[s_, l_, m_, n_, 0, opts___] := QNMRadial[s, l, m, n, 0.0, opts];


(* ::Section::Closed:: *)
(*QNMRadialFunction*)


SyntaxInformation[QNMRadialFunction] =
 {"ArgumentsPattern" -> {_}};


SetAttributes[QNMRadialFunction, {NHoldAll}];


(* ::Subsection::Closed:: *)
(*Output format*)


QNMRadialFunction /:
MakeBoxes[qnmrf: QNMRadialFunction[assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", assoc["s"]}], "  ",
                  BoxForm`SummaryItem[{"l: ", assoc["l"]}], "  ",
                  BoxForm`SummaryItem[{"m: ", assoc["m"]}], "  ",
                  BoxForm`SummaryItem[{"n: ", assoc["n"]}]}],
             BoxForm`SummaryItem[{"a: ", assoc["a"]}]};
  extended = {BoxForm`SummaryItem[{"Frequency: ", assoc["\[Omega]"]}],
              BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Coordinates: ", assoc["Coordinates"]}]
              };
  BoxForm`ArrangeSummaryBox[
    QNMRadialFunction,
    qnmrf,
    None,
    summary,
    extended,
    form
  ]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


QNMRadialFunction[assoc_][key_String] := assoc[key];


QNMRadialFunction[assoc_]["RadialFunction"] := Missing["KeyAbsent", "RadialFunction"];


Keys[m_QNMRadialFunction] ^:= DeleteCases[Join[Keys[m[[-1]]], {}], "RadialFunction"];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


outsideDomainQ[r_, rmin_, rmax_] := Min[r]<rmin || Max[r]>rmax;


QNMRadialFunction[assoc_][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[QNMRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[assoc["RadialFunction"][r], InterpolatingFunction::dmval]
 ];


Derivative[n_][QNMRadialFunction[assoc_]][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[QNMRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[Derivative[n][assoc["RadialFunction"]][r], InterpolatingFunction::dmval]
 ];


(*Derivative[n_Integer/;n>1][qnmrf:(QNMRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_])][r0:(_?NumericQ|{_?NumericQ..})] :=
 Module[{Rderivs, R, r, i, res},
  Rderivs = D[R[r_], {r_, i_}] :> D[0(*FIXME*), {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Collect[D[Derivative[i - 1][R][r], r] /. Rderivs, {R'[r], R[r]}, Simplify];, {i, 2, n}];
  res = Derivative[n][R][r] /. {
    R'[r] -> qnmrf'[r0],
    R[r] -> qnmrf[r0], r -> r0};
  Clear[Rderivs, i];
  Remove[R, r];
  res
];*)


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{QNMFrequency, QNMFrequencyKN, QNMRadial, QNMRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];
