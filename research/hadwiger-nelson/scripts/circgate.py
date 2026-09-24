"""Circular colourings through a real character: c(x) = floor(5 * frac(phi(x))), phi in Hom(M, R/Z).

Proper exactly when frac(phi(u)) lies in [1/5, 4/5] for every unit u.  Coset colourings are the
points phi = psi/5 of order 5; there every unit sits at 1/5, 2/5, 3/5 or 4/5, so the least slack
s = min_u dist(frac(phi(u)), {1/5, 4/5} outside) is 0 there, and by the Stiemke certificates they
are isolated.  Any phi with s > 0 gives a whole open set of characters, irrational ones among them:
quasi-periodic colourings of the module that are nowhere coset and split 5e pairs.

MILP (HiGHS via scipy):  maximise s  subject to  n_u + 1/5 + s <= C_u . phi <= n_u + 4/5 - s,
phi in [0, 1]^r, n_u integer, one u per +-pair.

usage: python3 circgate.py <module.json> [time limit s] [k]   (k colours, default 5)"""
import sys, json, os, time, itertools
from fractions import Fraction as Fr
from math import gcd
import numpy as np
import sympy
from scipy.optimize import milp, LinearConstraint, Bounds
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
t0 = time.time()
path = sys.argv[1]; TL = float(sys.argv[2]) if len(sys.argv) > 2 else 600; KC = int(sys.argv[3]) if len(sys.argv) > 3 else 5
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(p, q) for p, q in xy[0]]), F.element([Fr(p, q) for p, q in xy[1]]))
if "units" in d: U = [mk(xy) for xy in d["units"]]
else:
    P = [mk(xy) for xy in d["points"]]; g = build_graph(P); Ud = {}
    for i, j in g.edges():
        for w in (P[j] - P[i], P[i] - P[j]): Ud.setdefault((round(w.fx, 9), round(w.fy, 9)), w)
    U = list(Ud.values())
if F.gens[0] == 3:                              # close under the 60-degree turn when sqrt(3) is there
    cc = [Fr(0)] * F.dim; cc[1] = Fr(1); S3 = F.element(cc)
    turn = lambda p: Point(p.x * Fr(1, 2) - p.y * S3 * Fr(1, 2), p.x * S3 * Fr(1, 2) + p.y * Fr(1, 2))
    seen = {(round(u.fx, 7), round(u.fy, 7)) for u in U}; n0 = len(U)
    for u in list(U):
        w = u
        for _ in range(5):
            w = turn(w); k_ = (round(w.fx, 7), round(w.fy, 7))
            if k_ not in seen: seen.add(k_); U.append(w)
else: n0 = len(U)
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
den = 1
for v in Ev:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev[:n0])})
B = echelon(E); r = len(B)
C0 = [coords(B, tuple(int(Fr(x) * den) for x in v)) for v in Ev]
assert all(Fr(t).denominator == 1 for c in C0 for t in c)
C0 = [[int(t) for t in c] for c in C0]
src = open(os.path.join(HERE, "allunits.py")).read()
exec(src[src.index("def lll_exact"):src.index("Gr, Ur = lll_exact(G)")])
wt = list(F._prod) * 2
G0 = [[Fr(sum(w * p * q for w, p, q in zip(wt, B[i], B[j]))) for j in range(r)] for i in range(r)]
_, T = lll_exact(G0)
Ti = sympy.Matrix(T).inv()
CU = np.array([[int(x) for x in (sympy.Matrix([c]) * Ti)] for c in C0], dtype=np.int64)
# optional: also require the multiples k*u (MULT="1,5" asks for phi and 5*phi both circular colourings)
MULT = [int(x) for x in os.environ.get("MULT", "1").split(",")]
if MULT != [1]:
    CU = np.vstack([k * CU for k in MULT])
    print(f"  requiring the multiples {MULT} of every unit: {len(CU)} vectors", flush=True)
# one representative per +-pair
reps, seenr = [], set()
for i in range(len(CU)):
    k_ = tuple(CU[i]); nk_ = tuple(-CU[i])
    if nk_ in seenr or k_ in seenr: continue
    seenr.add(k_); reps.append(i)
C = CU[reps].astype(float); m = len(reps)
lo_n = np.floor(np.minimum(C, 0).sum(axis=1)) - 1; hi_n = np.ceil(np.maximum(C, 0).sum(axis=1)) + 1
print(f"{os.path.basename(path)}: {len(U)} units ({m} directions), rank {r}; max |coef| {np.abs(CU).max()}   [{time.time()-t0:.0f}s]", flush=True)
# variables: phi (r), n (m), s (1).  Maximise s.
nv = r + m + 1
cobj = np.zeros(nv); cobj[-1] = -1.0
A1 = np.zeros((m, nv)); A1[:, :r] = C; A1[np.arange(m), r + np.arange(m)] = -1.0; A1[:, -1] = -1.0   # C phi - n - s >= 1/K
A2 = np.zeros((m, nv)); A2[:, :r] = C; A2[np.arange(m), r + np.arange(m)] = -1.0; A2[:, -1] = 1.0    # C phi - n + s <= (K-1)/K
cons = [LinearConstraint(A1, lb=1.0 / KC, ub=np.inf), LinearConstraint(A2, lb=-np.inf, ub=(KC - 1.0) / KC)]
lb = np.concatenate([np.zeros(r), lo_n, [0.0]]); ub = np.concatenate([np.ones(r), hi_n, [0.5]])
integrality = np.concatenate([np.zeros(r), np.ones(m), [0]])
PAIRCK = os.environ.get("PAIRCK")
if PAIRCK:
    # the pair (A, B) of a checkpoint in the same module: require frac(phi(B - A)) in [DELTA, 1 - DELTA],
    # i.e. a circular colouring that splits the pair on a set of positive density
    dk = json.load(open(PAIRCK)); DELTA = float(os.environ.get("DELTA", "0.001"))
    PA = mk(dk["points"][dk["A"]]); PB = mk(dk["points"][dk["B"]])
    w = PB - PA; vw = tuple(w.x.c) + tuple(w.y.c)
    assert all((Fr(x) * den).denominator == 1 for x in vw), "B - A has a denominator outside den"
    cw = coords(B, tuple(int(Fr(x) * den) for x in vw))
    assert all(Fr(t).denominator == 1 for t in cw), "B - A is not in the module (check den)"
    cw = np.array([int(x) for x in (sympy.Matrix([[int(t) for t in cw]]) * Ti)], dtype=float)
    lo_m = np.floor(np.minimum(cw, 0).sum()) - 1; hi_m = np.ceil(np.maximum(cw, 0).sum()) + 1
    nv2 = nv + 1
    cobj = np.concatenate([cobj, [0.0]])
    cons = [LinearConstraint(np.hstack([c.A, np.zeros((c.A.shape[0], 1))]), c.lb, c.ub) for c in cons]
    row = np.zeros((1, nv2)); row[0, :r] = cw; row[0, -1] = -1.0
    cons.append(LinearConstraint(row, lb=DELTA, ub=1.0 - DELTA))
    lb = np.concatenate([lb, [lo_m]]); ub = np.concatenate([ub, [hi_m]]); integrality = np.concatenate([integrality, [1]])
    print(f"  pair constraint: frac(phi(B - A)) in [{DELTA}, {1 - DELTA}]; |B - A|^2 = {w.fx**2 + w.fy**2:.4f}", flush=True)
res = milp(cobj, constraints=cons, integrality=integrality, bounds=Bounds(lb, ub), options={"time_limit": TL, "disp": False})
print(f"  status {res.status}: {res.message}", flush=True)
if res.x is not None:
    phi = res.x[:r]; s = res.x[r + m]
    vals = (C @ phi) % 1.0
    print(f"  best slack s = {s:.6f}; phi = {np.round(phi, 6)}")
    print(f"  frac(phi(u)) range [{vals.min():.4f}, {vals.max():.4f}];  5*phi on the basis: {np.round(5 * phi, 4)}")
    if hasattr(res, "mip_dual_bound") and res.mip_dual_bound is not None: print(f"  dual bound on s: {-res.mip_dual_bound:.6f}")
print(f"  [{time.time()-t0:.0f}s]")
