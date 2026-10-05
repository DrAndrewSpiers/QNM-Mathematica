# QNM

Copyright (c) 2024-2025 Black Hole Perturbation Toolkit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17114757.svg)](https://doi.org/10.5281/zenodo.17114757)

A Mathematica package for computing quasinormal mode solutions of the Teukolsky equation in Schwarzschild and Kerr spacetime, and quasinormal modes of Kerr-Newman black holes.

Further details can be found on the [package homepage](https://bhptoolkit.org/QNM-Mathematica).

### Kerr-Newman black holes

Quasinormal modes of Kerr-Newman black holes with spin `a` and charge `Q` (in units of the mass) are computed from the separated Einstein-Maxwell perturbation equations of P. Hintz, "Perturbations of Kerr-Newman black holes: separation and mode stability", [arXiv:2609.33661](https://arxiv.org/abs/2609.33661) (2026):

```Mathematica
QNMFrequencyKN[s, l, m, n, a, Q]
QNMRadialKN[s, l, m, n, a, Q]
```

- `|s| = 2` gives the gravitational-led and `|s| = 1` the electromagnetic-led mode; they reduce to the Kerr modes with spin weight `s` as `Q -> 0`. The dipole `l = 1` exists only for the electromagnetic-led family, and `s = 0` is not supported.
- The frequency does not depend on the sign of `s` or of `Q`, nor on Hintz's spin system. For `Q = 0`, `QNMFrequencyKN` returns `QNMFrequency[-Abs[s], l, m, n, a]`.
- Small charges are supported; frequencies have been checked down to `Q = 10^-14`, where they agree with the Kerr values (the difference is of order `Q^2`). For small `Q` the gravitational-led radial function has its asymptotic form only for `r >> 1/Q^2` in the default spin system; use `Method -> {"HintzSeparated", "SpinSystem" -> -1}` there.
- `QNMRadialKN` returns Hintz's radial function `w` (or `w♯`), normalised to `w = 1` at the outer horizon, as a `QNMRadialFunction`.

```Mathematica
In[1]:= QNMFrequencyKN[-2, 2, 2, 0, 0.6, 0.5]
Out[1]= 0.536127 - 0.0811784 I

In[2]:= w = QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5]; w[10]
Out[2]= -6.8711 - 13.6778 I
```

### License

This code is distributed under the MIT License. Details can be found in the LICENSE file.
