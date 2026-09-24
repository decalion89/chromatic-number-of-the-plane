"""Sample the cell of F(M1) around a known point and test every sample's extensions to M'.

The cell: fix the integer parts n_u of phi(u) at the known point; the cell is the polytope
{phi : n_u + 1/5 <= C_u . phi <= n_u + 4/5}.  Vertices come from LPs with random objectives;
samples are random convex combinations.  Each sample is pushed through circextend's enumeration of
the |det J| extensions to M' (vectorised, float, exact re-check of any hit).

usage: python3 cellsample.py <M1.json> <Mprime.json> <samples> <phi_1..phi_r>"""
import sys, os, time, json
from fractions import Fraction as Fr
from math import floor
import numpy as np, sympy
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
M1, MP, NS = sys.argv[1], sys.argv[2], int(sys.argv[3]); phi0 = [Fr(x) for x in sys.argv[4:]]
# reuse circextend's set-up (bases, J, Smith form, extension grid) by running its top part
src = open(os.path.join(HERE, "circextend.py")).read()
top = src[:src.index("Lphi = [sum(L[i, j] * phi[j] for j in range(r)) for i in range(r)]")]
sys.argv = ["circextend.py", M1, MP] + [str(x) for x in phi0]
exec(top)
CU1 = g1["CU"].astype(float); r = g1["r"]
C1 = g1["CU"]
nfix = [floor(sum(Fr(int(c)) * p for c, p in zip(cu, phi0))) for cu in C1]
A_ub = np.vstack([C1.astype(float), -C1.astype(float)])
b_ub = np.concatenate([np.array([n + 0.8 for n in nfix]), -np.array([n + 0.2 for n in nfix])])
rng = np.random.default_rng(1)
verts = []
for s_ in range(60):
    res = linprog(rng.normal(size=r), A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * r, method="highs")
    if res.status == 0: verts.append(res.x)
verts = np.unique(np.round(np.array(verts), 9), axis=0)
print(f"cell vertices found: {len(verts)}; spread per coordinate: {np.round(verts.max(0) - verts.min(0), 4)}   [{time.time()-t0:.0f}s]", flush=True)
Rn = np.array([[float(R[i, j]) for j in range(r)] for i in range(r)])
CUp = gp["CU"].astype(float)
Lf = np.array([[float(L[i, j]) for j in range(r)] for i in range(r)])
dgf = np.array([float(x) for x in dg])
grid = np.stack([g_.ravel() for g_ in np.meshgrid(*[np.arange(abs(int(x)), dtype=float) for x in dg], indexing="ij")], axis=1)
best_all = -9
for s_ in range(NS):
    w = rng.dirichlet(np.ones(len(verts)) * 0.5)
    ph = w @ verts
    # extensions: z = (L phi + j) / d, y = R z
    Z = (Lf @ ph + grid) / dgf
    Y = Z @ Rn.T
    vals = (Y @ CUp.T) % 1.0
    sl = np.minimum(vals - 0.2, 0.8 - vals).min(axis=1)
    b = sl.max(); best_all = max(best_all, b)
    if b >= -1e-9:
        print(f"  sample {s_}: an extension works (slack {b:.6f})!  phi = {np.round(ph, 6)}", flush=True)
        break
print(f"{NS} samples; best extension slack over all samples: {best_all:.6f}   [{time.time()-t0:.0f}s]")
