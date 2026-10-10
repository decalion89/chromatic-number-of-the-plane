"""Independent second check of a clique-LP certificate written by certify_lp2.py (NAME.cert.json).

Different code path: sympy for exact clique distances and radii, mpmath (60 digits) for Bessel values, a fixed
Lipschitz constant computed from scratch, an adaptive sweep of its own, and the tail bound |J0(x)| <= sqrt(2/(pi x)).
This is a high-precision floating-point check, not interval arithmetic; check3.py is the rigorous independent check.
usage: python3 check2.py NAME.cert.json k
"""
import sys, json, itertools
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
from shapes import SHAPES
mp.mp.dps = 60

def S(e):
    return sp.sympify(e, locals={'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'pi': sp.pi})

C = json.load(open(sys.argv[1])); k = int(sys.argv[2])
D = [S(e) for e in C['D'].split(';')]
z0 = Fr(C['z0']); z = [Fr(v) for v in C['z']]
terms = [(z[i], D[i]) for i in range(len(D))]        # (coefficient, radius) pairs, radius exact sympy
ysum = Fr(0); m = Fr(0)
for c in C['cliques']:
    W = [(S(x), S(y)) for x, y in SHAPES[c['shape']]]
    for a, b in itertools.combinations(range(len(W)), 2):
        d2 = (W[a][0] - W[b][0])**2 + (W[a][1] - W[b][1])**2
        assert any(sp.simplify(sp.expand(d2 - dd**2)) == 0 for dd in D), 'clique'   # symbolic, no numerics
    cx = sum(w[0] for w in W) / len(W); cy = sum(w[1] for w in W) / len(W)
    px, py = cx + sp.Rational(c['ox']), cy + sp.Rational(c['oy'])
    y = Fr(c['y']); assert y >= 0
    for wx, wy in W:
        terms.append((y, sp.sqrt((wx - px)**2 + (wy - py)**2)))
    ysum += y; m += y * len(W)
R = [mp.mpf(str(sp.N(r, 70))) for _, r in terms]
Cf = [mp.mpf(cf.numerator) / cf.denominator for cf, _ in terms]
Z0 = mp.mpf(z0.numerator) / z0.denominator
L = sum(abs(c) * r for c, r in zip(Cf, R)) / mp.sqrt(2) * (1 + mp.mpf('1e-30'))
g = lambda t: Z0 + sum(c * mp.besselj(0, r * t) for c, r in zip(Cf, R))
T = mp.mpf(C['T'])
t = mp.mpf(0); n = 0; gmin = None
while t < T:
    v = g(t)
    if v <= 0:
        raise SystemExit(f'FAIL: g({t}) = {v}')
    gmin = v if gmin is None or v < gmin else gmin
    t += v / L * mp.mpf('0.99'); n += 1
tail = sum(abs(c) * mp.sqrt(2 / (mp.pi * r * T)) for c, r in zip(Cf, R))
assert Z0 - tail > 0, 'tail'
s = z0 + sum(z) + m
B = (z0 + ysum) / s
print(f"check2: {n} steps on [0, {C['T']}], min g {mp.nstr(gmin, 6)}, tail {mp.nstr(tail, 6)} < z0 {mp.nstr(Z0, 6)}")
print(f"check2: B = {float(B):.8f}, 1/B = {float(1/B):.5f}, B < 1/{k}: {B < Fr(1, k)}")
