"""Referee E: exact check of Lemma 'radii' (paper, Section 9) which Lemma 12 uses for the bound R.
(i) P_1 (s=1) = {x : conj(x) rho^j in B_1, |j|<=1}: vertices and max modulus^2 (claim: 12-gon, max |x|^2 = 10/9).
(ii) X((1+i)/3) at level 5 = {y : eta_{j,e} (conj(y) rho^j)_e <= 1, |j|<=1, e=1,2} with eta_{-1}=1+i, eta_0=1-i,
eta_1=-1-i (recomputed here from Lemma 'values'(iii)); claim: hexagon with vertices -5/4, -5i/4, -5/7(1+i),
1-i/2, 1+i, -1/2+i, max modulus sqrt 2.  Also recheck: for theta in (2/7,1/3], R = max(s sqrt10/3, sqrt2 delta) < 0.23."""
from fractions import Fraction as Fr
import itertools
from math import floor

def conj_rho(y, j):
    re, im = y[0], -y[1]
    r = (Fr(3, 5), Fr(4, 5)) if j >= 0 else (Fr(3, 5), Fr(-4, 5))
    for _ in range(abs(j)):
        re, im = re * r[0] - im * r[1], re * r[1] + im * r[0]
    return (re, im)

def lin(j, e):
    """coefficients (a,b) with (conj(y) rho^j)_e = a*y1 + b*y2"""
    u = conj_rho((Fr(1), Fr(0)), j)
    v = conj_rho((Fr(0), Fr(1)), j)
    return (u[e], v[e])

def vertices(cons):
    """cons: list of (a,b,c) meaning a*y1+b*y2 <= c. Return vertices (exact) of the polygon."""
    vs = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(cons, 2):
        d = a1 * b2 - a2 * b1
        if d == 0:
            continue
        x = (c1 * b2 - c2 * b1) / d
        y = (a1 * c2 - a2 * c1) / d
        if all(a * x + b * y <= c for (a, b, c) in cons):
            vs.add((x, y))
    return sorted(vs)

# (i)
cons = []
for j in (-1, 0, 1):
    for e in (0, 1):
        a, b = lin(j, e)
        cons.append((a, b, Fr(1)))
        cons.append((-a, -b, Fr(1)))
V = vertices(cons)
print("(i) P_1 for s=1:", len(V), "vertices:", [(str(x), str(y)) for x, y in V])
print("    max |x|^2 =", max(x * x + y * y for x, y in V), "(claim 10/9)")
# (ii)
N = 5
h = (Fr(1, 2), Fr(1, 2))
eps = (Fr(1, 3), Fr(1, 3))
etas = {}
for j in (-1, 0, 1):
    w = conj_rho((N * eps[0], N * eps[1]), j)
    w = (w[0] - floor(w[0]), w[1] - floor(w[1]))
    etas[j] = (6 * (w[0] - h[0]), 6 * (w[1] - h[1]))
print("(ii) eta_j((1+i)/3) at N=5:", {j: (str(a), str(b)) for j, (a, b) in etas.items()})
cons = []
for j in (-1, 0, 1):
    for e in (0, 1):
        a, b = lin(j, e)
        cons.append((etas[j][e] * a, etas[j][e] * b, Fr(1)))
V = vertices(cons)
print("    X vertices:", [(str(x), str(y)) for x, y in V])
print("    max |y|^2 =", max(x * x + y * y for x, y in V), "(claim 2)")
# R bound: s < 3/14, delta < 1/21 for theta > 2/7; R^2 < max(10/9 * 9/196, 2/441)
print("R^2 sup over theta>2/7:", max(Fr(10, 9) * Fr(9, 196), Fr(2, 441)), "= %.6f ; 0.23^2 = 0.0529" % float(max(Fr(10, 9) * Fr(9, 196), Fr(2, 441))))
