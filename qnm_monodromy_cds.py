# -*- coding: utf-8 -*-
"""
Asymptotic gravitational quasinormal modes of the C-de Sitter (CdS) sector of the
vacuum Plebanski-Demianski family, from the monodromy quantization condition

    F(omega) = cosh(pi*omega/k_+ - pi*omega/k_c) + 3*cosh(pi*omega/k_+ + pi*omega/k_c) = 0 .

Time convention e^{-i omega t}: damped modes have Im(omega) < 0.
The file is organized in # %% cells (Spyder / VS Code): run CONFIGURATION and SETUP once,
then SOLVE, TABLE and PLOT independently.
"""

# %% CONFIGURATION ------------------------------------------------------------------
m      = '1'        # mass parameter
alpha  = '0.1'      # acceleration
Lam    = '0.01'     # cosmological constant (Lambda >= 0)
V_MAX  = 1.75       # modes are searched for 0 < -Im(omega) < V_MAX
DPS    = 50         # decimal digits of precision (mpmath)
CSV    = 'cds_qnm_modes.csv'
FIG    = 'modes_CdS.png'

# %% SETUP: horizons and surface gravities ---------------------------------------------
import csv
import mpmath as mp
mp.mp.dps = DPS
m, alpha, Lam = mp.mpf(m), mp.mpf(alpha), mp.mpf(Lam)

def Q(r):  return (r**2 - 2*m*r)*(1 - alpha**2*r**2) - Lam*r**4/3
def dQ(r): return mp.diff(Q, r)

# Q(r)/r is a cubic: its three roots are r_n < 0 < r_+ < r_c
coeffs = [-(alpha**2 + Lam/3), 2*m*alpha**2, 1, -2*m]
raw = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
if any(abs(mp.im(z)) > mp.mpf(10)**(-DPS//2) for z in raw):
    raise ValueError('Q(r)/r has complex roots: the horizon structure assumed here does not exist.')
r_n, r_p, r_c = sorted(mp.re(z) for z in raw)
if not (r_n < 0 < r_p < r_c):
    raise ValueError('Expected r_n < 0 < r_+ < r_c (no extremal or Nariai configuration).')

kappa = {name: dQ(r)/(2*r**2) for name, r in (('n', r_n), ('+', r_p), ('c', r_c))}
k_n, k_p, k_c = kappa['n'], kappa['+'], kappa['c']
lam1 = mp.pi*abs(1/k_p - 1/k_c)
lam2 = mp.pi*abs(1/k_p + 1/k_c)
u_star = mp.log(3)/(lam1 - lam2)

print(f'r_n = {mp.nstr(r_n, 10)}   r_+ = {mp.nstr(r_p, 10)}   r_c = {mp.nstr(r_c, 10)}')
print(f'k_n = {mp.nstr(k_n, 10)}   k_+ = {mp.nstr(k_p, 10)}   k_c = {mp.nstr(k_c, 10)}')
print(f'check 1/k_n + 1/k_+ + 1/k_c = {mp.nstr(1/k_n + 1/k_p + 1/k_c, 3)}')
print(f'u* = ln3 max(k_+,|k_c|)/2pi = {mp.nstr(u_star, 10)}')

def F(w):
    return mp.cosh(mp.pi*w/k_p - mp.pi*w/k_c) + 3*mp.cosh(mp.pi*w/k_p + mp.pi*w/k_c)

# %% SOLVE: zeros of F in the Polya strips ---------------------------------------------
# (1) imaginary axis: F(iv) = cos(lam1 v) + 3 cos(lam2 v) is real -> sign changes + bisection
def G(v): return mp.cos(lam1*v) + 3*mp.cos(lam2*v)
step = (mp.pi/lam1)/200
axis, v = [], mp.mpf(0)
while v > -V_MAX:
    v2 = v - step
    if G(v)*G(v2) < 0:
        axis.append(mp.mpc(0, mp.findroot(G, (v2, v), solver='bisect', tol=mp.mpf(10)**(-DPS+10))))
    v = v2

# (2) strip at u = +u*: start from the zeros of the two dominant terms
osc, k = [], 0
while True:
    w0 = (mp.log(3) - 1j*mp.pi*(2*k + 1))/(lam1 - lam2)
    if -mp.im(w0) > V_MAX + 1: break
    w = mp.findroot(F, w0, tol=mp.mpf(10)**(-DPS+10), maxsteps=100)
    if -V_MAX < mp.im(w) < 0 and all(abs(w - z) > 1e-12 for z in osc):
        osc.append(w)
    k += 1

modes = sorted([('I', z) for z in axis] + [('O', z) for z in osc], key=lambda t: -mp.im(t[1]))
print(f'{len(axis)} purely imaginary + {len(osc)} oscillating modes with -Im(omega) < {V_MAX}')
print(f'expected from the average spacings: ~{mp.nstr(V_MAX/k_n, 3)} and ~{mp.nstr(V_MAX/max(k_p, -k_c), 3)}')

# %% TABLE -----------------------------------------------------------------------------
with open(CSV, 'w', newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['n', 'family', 'Re_omega', 'Im_omega', 'abs_F'])
    print(f"{'n':>3} {'':>2} {'Re(omega)':>12} {'Im(omega)':>12} {'|F|':>9}")
    for n, (fam, z) in enumerate(modes):
        re_, im_, res = mp.nstr(mp.re(z), 10), mp.nstr(mp.im(z), 10), mp.nstr(abs(F(z)), 2)
        wr.writerow([n, fam, re_, im_, res])
        print(f'{n:>3} {fam:>2} {re_:>12} {im_:>12} {res:>9}')
print(f'saved {CSV}')

# %% PLOT ------------------------------------------------------------------------------
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(5.2, 4.4))
xo = [float(mp.re(z)) for f, z in modes if f == 'O']; yo = [float(mp.im(z)) for f, z in modes if f == 'O']
yi = [float(mp.im(z)) for f, z in modes if f == 'I']
ax.plot(xo, yo, '-o', color='#1f77b4', ms=4, lw=1, label='oscillating')
ax.plot([0]*len(yi), yi, '^', color='#ff7f0e', ms=5, label='purely imaginary')
ax.axvline(float(u_star), color='0.4', ls=':', lw=0.9, label=r'$u^*$')
ax.set_xlabel(r'$\mathrm{Re}(m\omega)$'); ax.set_ylabel(r'$\mathrm{Im}(m\omega)$')
ax.grid(alpha=0.25, ls=':'); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(FIG, dpi=300)
print(f'saved {FIG}')
plt.show()
