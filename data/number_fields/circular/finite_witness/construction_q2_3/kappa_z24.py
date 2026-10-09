"""kappa_z24.py SPEC...: kappa (float MILP) of U = mu_24 * {g^l : |l| <= L} for rotations g in Q(zeta_24) given by
name: om7 = (8 + 5 omega)/7, s2 = (1 + 2 sqrt-2)/3, s6 = (5 + 2 sqrt-6)/7, i5 = (3 + 4i)/5; SPEC is name:L."""
import sys
from fractions import Fraction as Fr
from math import lcm
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from sympy import Matrix
N = 24; D = 8
PHI = [1, 0, 0, 0, -1, 0, 0, 0, 1]          # zeta^8 - zeta^4 + 1


def mul(a, b):
    r = [Fr(0)] * (2 * D - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    for k in range(2 * D - 2, D - 1, -1):       # reduce zeta^k = zeta^(k-8) * zeta^8 = zeta^(k-8) (zeta^4 - 1)
        c = r[k]
        if c:
            r[k] = 0; r[k - 4] += c; r[k - 8] -= c
    return r[:D]


def zpow(j):
    e = [Fr(0)] * D; e[0] = Fr(1)
    z = [Fr(0)] * D; z[1] = Fr(1)
    for _ in range(j % N):
        e = mul(e, z)
    return e


def inv_unit(g):
    """inverse of a norm-1 element = complex conjugate: zeta -> zeta^-1"""
    r = [Fr(0)] * D
    for i, c in enumerate(g):
        if c:
            t = zpow(-i)
            r = [x + c * y for x, y in zip(r, t)]
    return r


z = lambda j: zpow(j)
omega = z(8); i_ = z(6); s3 = [a + b for a, b in zip(z(2), z(-2))]       # zeta^2 + zeta^-2 = 2 cos 30 = sqrt3
sqrtm3 = [a - b for a, b in zip(omega, inv_unit(omega))]                  # omega - omega^-1 = i sqrt3
sqrt2 = [a + b for a, b in zip(z(3), z(-3))]                              # 2 cos 45 = sqrt2
sqrtm2 = mul(i_, sqrt2); sqrtm6 = mul(i_, mul(sqrt2, s3))
lin = lambda *terms: [sum(c * v[k] for c, v in terms) for k in range(D)]
ONE = z(0)
G = {"om7": lin((Fr(8, 7), ONE), (Fr(5, 7), omega)),
     "s2": lin((Fr(1, 3), ONE), (Fr(2, 3), sqrtm2)),
     "s6": lin((Fr(5, 7), ONE), (Fr(2, 7), sqrtm6)),
     "i5": lin((Fr(3, 5), ONE), (Fr(4, 5), i_))}
for k, g in G.items():
    assert mul(g, inv_unit(g)) == ONE, k
specs = [(s.split(":")[0], int(s.split(":")[1])) for s in sys.argv[1:]]
gens = [ONE]
for name, L in specs:
    g = G[name]; gi = inv_unit(g)
    new = []
    for v in gens:
        p = v
        for l in range(1, L + 1):
            p = mul(p, g); new.append(p)
        p = v
        for l in range(1, L + 1):
            p = mul(p, gi); new.append(p)
    gens = gens + new
U = []
seen = set()
for v in gens:
    for j in range(N // 2):
        w = tuple(mul(v, z(j)))
        if w not in seen and tuple(-x for x in w) not in seen:
            seen.add(w); U.append(w)
den = lcm(*[x.denominator for u in U for x in u])
Mi = Matrix([[int(x * den) for x in u] for u in U])        # rows: vectors (integers)
# Z-basis of the row lattice by HNF: coordinates of each vector in that basis
H = Mi.T.echelon_form()   # placeholder to keep sympy imported
from sympy.matrices.normalforms import hermite_normal_form
Hn = hermite_normal_form(Mi.T)                              # columns span the lattice
B = np.array(Hn.T.tolist(), dtype=object)                   # rows = basis vectors
Bf = Matrix(B.tolist())
coords = [list(Bf.T.solve(Matrix([int(x * den) for x in u]))) for u in U]
C = np.array([[int(c) for c in row] for row in coords], dtype=float)
n = len(U); r = C.shape[1]
nv = r + n + 1
cobj = np.zeros(nv); cobj[-1] = -1
A = []; lo = []; hi = []
for idx in range(n):
    row = np.zeros(nv); row[:r] = C[idx]; row[r + idx] = -1
    r1 = row.copy(); r1[-1] = -1; A.append(r1); lo.append(0); hi.append(np.inf)
    r2 = row.copy(); r2[-1] = 1; A.append(r2); lo.append(-np.inf); hi.append(1)
kmax = int(np.abs(C).sum(axis=1).max()) + 2
bounds = Bounds([0] * r + [-kmax] * n + [0], [1] * r + [kmax] * n + [0.5])
res = milp(cobj, constraints=LinearConstraint(np.array(A), lo, hi), integrality=np.array([0] * r + [1] * n + [0]),
           bounds=bounds, options={"time_limit": 300})
print(f"{sys.argv[1:]}: {n} vectors (up to sign), rank {r}, kappa ~ {res.x[-1]:.6f}" if res.x is not None else res.message)
