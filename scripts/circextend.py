"""Does a circular colouring of a module M1 extend to a bigger module M' (e.g. its lambda-closure)?

phi is given on circgate.py's LLL basis of M1.  The extensions phi' of phi to M' are the solutions of
J y = phi (mod Z^r), J the inclusion matrix of M1's basis in M''s basis: |det J| of them, enumerated
through the Smith form of J.  Each is tested on every unit of M' (circgate.py's basis of M').

usage: python3 circextend.py <M1.json> <Mprime.json> <phi_1> ... <phi_r>"""
import sys, json, os, time, itertools
from fractions import Fraction as Fr
from math import gcd, floor
import numpy as np, sympy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
t0 = time.time()
M1, MP = sys.argv[1], sys.argv[2]; phi = [Fr(x) for x in sys.argv[3:]]
src = open(os.path.join(HERE, "circgate.py")).read()
pre = src[:src.index("# one representative per +-pair")]
def basis_of(path):
    g = {"__file__": os.path.join(HERE, "circgate.py")}
    sys.argv = ["circgate.py", path, "1"]
    exec(pre, g)
    return g
g1 = basis_of(M1); gp = basis_of(MP)
F = gp["F"]; r = gp["r"]
assert tuple(g1["F"].gens) == tuple(F.gens) and g1["r"] == r
# LLL basis vectors of M1 as field vectors: rows of T1 * B1 / den1
def lll_rows(g):
    T = sympy.Matrix(g["T"]); Bm = sympy.Matrix([[int(x) for x in b] for b in g["B"]])
    return (T * Bm) / g["den"]
R1 = lll_rows(g1); Rp = lll_rows(gp)
# coordinates of R1's rows in Rp's basis: J = R1 * Rp^+ (exact via pivot columns)
cols = []
for j in range(Rp.shape[1]):
    if sympy.Matrix.hstack(*([Rp[:, k] for k in cols] + [Rp[:, j]])).rank() > len(cols): cols.append(j)
    if len(cols) == r: break
J = R1[:, cols] * Rp[:, cols].inv()
assert J * Rp == R1, "M1 is not inside M'"
assert all(x.q == 1 for x in J), "M1's basis is not integral in M'"
detJ = abs(J.det())
print(f"[M' : M1] = {detJ}   [{time.time()-t0:.0f}s]", flush=True)
# extensions: y with J y = phi + k  (column convention: phi(b1_i) = sum_j J_ij y_j)
from sympy.matrices.normalforms import smith_normal_form
# Smith form with transforms by elimination (small 8x8)
def snf(A):
    A = sympy.Matrix(A); m, n_ = A.shape; L = sympy.eye(m); R = sympy.eye(n_); k = 0
    while k < min(m, n_):
        best = None
        for i in range(k, m):
            for j in range(k, n_):
                if A[i, j] != 0 and (best is None or abs(A[i, j]) < abs(A[best[0], best[1]])): best = (i, j)
        if best is None: break
        i, j = best; A.row_swap(k, i); L.row_swap(k, i); A.col_swap(k, j); R.col_swap(k, j)
        done = True
        for i in range(k + 1, m):
            q = A[i, k] // A[k, k]
            if q: A[i, :] = A[i, :] - q * A[k, :]; L[i, :] = L[i, :] - q * L[k, :]
            if A[i, k] != 0: done = False
        for j in range(k + 1, n_):
            q = A[k, j] // A[k, k]
            if q: A[:, j] = A[:, j] - q * A[:, k]; R[:, j] = R[:, j] - q * R[:, k]
            if A[k, j] != 0: done = False
        if done: k += 1
    return A, L, R
D, L, R = snf(J)                      # L J R = D  ->  J = L^-1 D R^-1
dg = [int(D[i, i]) for i in range(r)]
print(f"  Smith invariants of J: {dg}", flush=True)
# J y = phi + k  <=>  D (R^-1 y) = L (phi + k).  Let z = R^-1 y: z_i = (L(phi + k))_i / d_i.
# Distinct y mod Z^r: z_i ranges over (L phi)_i/d_i + j/d_i, j = 0..|d_i|-1 (k absorbed).
Lphi = [sum(L[i, j] * phi[j] for j in range(r)) for i in range(r)]
Rn = np.array([[float(R[i, j]) for j in range(r)] for i in range(r)])
CUp = gp["CU"].astype(float)                          # units of M' in M''s LLL coordinates
choices = [[(Fr(Lphi[i]) + j) / dg[i] for j in range(abs(dg[i]))] for i in range(r)]
tot = int(np.prod([len(c) for c in choices]))
print(f"  {tot} extensions to test on {len(CUp)} units   [{time.time()-t0:.0f}s]", flush=True)
# vectorised over the product of choices, in float (then exact re-check of any hit)
best = (-1, None)
grids = np.meshgrid(*[np.array([float(x) for x in c]) for c in choices], indexing="ij")
Z = np.stack([g_.ravel() for g_ in grids], axis=1)            # tot x r
Y = Z @ Rn.T                                                   # y = R z
vals = (Y @ CUp.T) % 1.0                                       # tot x units
sl = np.minimum(vals - 0.2, 0.8 - vals).min(axis=1)
i = int(np.argmax(sl))
print(f"  best extension slack (float): {sl[i]:.6f};  extensions with slack >= -1e-9: {(sl >= -1e-9).sum()}   [{time.time()-t0:.0f}s]", flush=True)
if sl[i] >= -1e-9:
    z = [choices[k][int(np.argmin(np.abs(np.array([float(x) for x in choices[k]]) - Z[i, k])))] for k in range(r)]
    y = [sum(Fr(int(R[a, b])) * z[b] for b in range(r)) for a in range(r)]
    ex = min(min((sum(Fr(int(c)) * yy for c, yy in zip(cu, y)) % 1) - Fr(1, 5), Fr(4, 5) - (sum(Fr(int(c)) * yy for c, yy in zip(cu, y)) % 1)) for cu in gp["CU"])
    print(f"  exact slack of that extension: {ex}   phi' = {[str(v) for v in y]}")
