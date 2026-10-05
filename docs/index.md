{% include head.html %}

<p>
 <h1 style="display:inline">QNM</h1> <span style="float:right;"><a href="https://bhptoolkit.org/mathematica-install.html" class = "code_btn">Install this package!</a></span>
</p>

A Mathematica package for computing quasinormal mode solutions of the Teukolsky equation in Schwarzschild and Kerr spacetime. The package introduces two new functions:
```Mathematica
QNMFrequency[s, l, m, n, a]
QNMRadial[s, l, m, n, a]
```
where  
$s$ - is the spin-weight of the field  
$l$ - is the multipolar index  
$m$ - is the azimuthal index  
$n$ - is the overtone number  
$a$ - is the black hole spin parameter

### Kerr-Newman black holes

Quasinormal modes of Kerr-Newman black holes are computed from the separated Einstein-Maxwell perturbation equations of [P. Hintz, "Perturbations of Kerr-Newman black holes: separation and mode stability", arXiv:2609.33661 (2026)](https://arxiv.org/abs/2609.33661). Two further functions are provided:
```Mathematica
QNMFrequencyKN[s, l, m, n, a, Q]
QNMRadialKN[s, l, m, n, a, Q]
```
where $Q$ is the black hole charge (in units of the mass, with $a^2 + Q^2 < 1$). The spin weight selects the family of coupled gravito-electromagnetic modes: $|s|=2$ gives the gravitational-led and $|s|=1$ the electromagnetic-led mode, which reduce to the corresponding Kerr modes as $Q \to 0$. The dipole $l=1$ exists only for the electromagnetic-led family, and $s=0$ is not supported. The frequency does not depend on the sign of $s$ or of $Q$, nor on Hintz's spin system, and for $Q=0$ `QNMFrequencyKN` returns `QNMFrequency[-Abs[s], l, m, n, a]`. Small charges are supported: frequencies have been checked down to $Q = 10^{-14}$, where they agree with the Kerr values up to corrections of order $Q^2$. `QNMRadialKN` returns Hintz's radial function $w$ (normalised to 1 at the outer horizon) as a `QNMRadialFunction`.

```Mathematica
In[1]:= QNMFrequencyKN[-2, 2, 2, 0, 0.6, 0.5]
Out[1]= 0.536127 - 0.0811784 I

In[2]:= w = QNMRadialKN[-2, 2, 2, 0, 0.6, 0.5]; w[10]
Out[2]= -6.8711 - 13.6778 I
```

### Further examples

Additional documentation can be found in the Mathematica documentation included with the package.

## Citing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17114757.svg)](https://doi.org/10.5281/zenodo.17114757)

In addition to acknowledging the Black Hole Perturbation Toolkit as suggested on the [front page](https://bhptoolkit.org) we also recommend citing the specific package version you use via the citation information on the package’s Zenodo page linked from the above DOI.

## Authors and contributors

Barry Wardell, Christiana Pantelidou, Brad Cownden, Jake Mac Uilliam, Conor O'Toole, Rodrigo Macedo, Jamil Assaad
