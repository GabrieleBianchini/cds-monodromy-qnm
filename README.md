# cds-monodromy-qnm

### Asymptotic Gravitational Quasinormal Modes of the C-de Sitter Black Hole

This project computes the highly damped quasinormal mode (QNM) spectrum of gravitational perturbations ($s=+2$) of a uniformly accelerated black hole in a de Sitter background, the C-de Sitter (CdS) sector of the vacuum Plebański–Demiański family, from the analytic quantization condition obtained with the monodromy method. All computations use arbitrary-precision arithmetic (`mpmath`).

The calculation is part of:

- G. Bianchini, S. Gambino, O. Luongo, *"Plebański-Demiański metric in vacuum: a study of gravitational quasinormal modes"* (in preparation).

The quantization condition generalizes the Schwarzschild-de Sitter result of:

- V. Cardoso, J. Natário, R. Schiappa, *"Asymptotic quasinormal frequencies for black holes in non-asymptotically flat spacetimes"*, J. Math. Phys. **45**, 4698 (2004). [arXiv:hep-th/0403132](https://arxiv.org/abs/hep-th/0403132)

which in turn extends the Schwarzschild computation of:

- L. Motl, A. Neitzke, *"Asymptotic black hole quasinormal frequencies"*, Adv. Theor. Math. Phys. **7**, 307 (2003). [arXiv:hep-th/0301173](https://arxiv.org/abs/hep-th/0301173)

The metric function of the CdS sector is

```
Q(r) = (r^2 - 2 m r)(1 - alpha^2 r^2) - (Lambda/3) r^4
```

with mass `m`, acceleration `alpha` and cosmological constant `Lambda >= 0`. Besides `r = 0`, `Q(r)` has three real roots `r_n < 0 < r_+ < r_c`: a negative root, the event horizon and the cosmological-acceleration horizon. For `Lambda = 0` the metric reduces to the C-metric, with `r_n = -1/alpha`, `r_+ = 2m`, `r_c = 1/alpha`.

---

## Code Description

### Surface gravities

The tortoise coordinate is `dr_*/dr = r^2/Q(r)`, and the surface gravity of each root is

```python
kappa = {name: dQ(r)/(2*r**2) for name, r in (('n', r_n), ('+', r_p), ('c', r_c))}
```

with `kappa_c < 0`. Since `r^2/Q(r)` has no pole at the origin and no residue at infinity, the three surface gravities satisfy `1/k_n + 1/k_+ + 1/k_c = 0`, which the code prints as a check.

### Quantization condition

The master equation, written in the Teukolsky form for the Weyl scalar of spin weight 2, contains a term linear in `omega`. Near `r = 0` this term turns the local equation into a Coulomb wave equation, whose solutions connect as those of the Bessel equation with index `j = 2`. The monodromy around the contour enclosing `r_+` and `r_c` then gives

```python
def F(w):
    return mp.cosh(mp.pi*w/k_p - mp.pi*w/k_c) + 3*mp.cosh(mp.pi*w/k_p + mp.pi*w/k_c)
```

with the time convention `e^{-i omega t}`, so that damped modes have `Im(omega) < 0`. Since `F` is even, modes with `Re(omega) != 0` come in pairs `omega`, `-conj(omega)`; only those with `Re(omega) >= 0` are reported.

### Locating all the roots

A local root finder alone does not guarantee that no root is missed. `F` is an exponential polynomial, and by Pólya's theorem its zeros lie in vertical strips of the complex plane fixed by its coefficients and exponents. Writing `lam1 = pi|1/k_+ - 1/k_c|` and `lam2 = pi|1/k_+ + 1/k_c|`, the strips are:

- the imaginary axis, where `F(iv) = cos(lam1 v) + 3 cos(lam2 v)` is real. Its zeros are found exactly, by sign changes on a fine grid followed by bisection, and their average spacing is `pi/lam2 = k_n`;
- the lines `Re(omega) = +-u*`, with `u* = ln3/(lam1 - lam2) = ln3 max(k_+, |k_c|)/(2 pi)`. Here each root is refined with `mp.findroot`, starting from the zeros of the two dominant exponentials, `omega_k = (ln3 - i pi (2k+1))/(lam1 - lam2)`. Their average spacing is `max(k_+, |k_c|)`, and their real parts oscillate around `u*`.

The code compares the number of roots found in each strip with the number expected from these average spacings, and reports `|F|` at every root.

---

## Usage

```bash
pip install -r requirements.txt
```

Edit the `CONFIGURATION` cell at the top of `qnm_monodromy_cds.py` to set `m`, `alpha`, `Lambda` and the search window `V_MAX`, then run the `SETUP` cell and the `SOLVE`, `TABLE` and `PLOT` cells. The file is organized in `# %%` cells (Spyder / VS Code style), so each step can be re-run independently. The modes are written to `cds_qnm_modes.csv` and plotted in `modes_CdS.png`.

With the default values (`m = 1`, `alpha = 0.1`, `Lambda = 0.01`) the code reproduces the table of asymptotic frequencies reported in the paper.

---

## Known limitations

- The quantization condition holds only in the asymptotic regime `|omega| -> infinity`; the low overtones must be computed with other methods (e.g. Leaver's continued fractions).
- The code assumes three distinct real nonvanishing roots of `Q(r)` and stops otherwise. The extremal and near-Nariai limits (`r_+ ~ r_c`) are not covered.
- Rotating and NUT-charged members of the Plebański–Demiański family are not covered: there the Fourier number conjugate to the canonical Killing coordinate grows with `omega`, and the asymptotic problem acquires the structure of the Kerr one, which requires a separate analysis.
