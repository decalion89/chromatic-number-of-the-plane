"""task2_balls.py -- 4-colourability of B_k(D) u conj(B_k(D)) and of B_k(U), all U-edges.
usage: python3 task2_balls.py NAME [tlimit]
"""
import sys, time, json
import numpy as np
from flat import *
from satutil import *
from pysat.solvers import Solver

D, CD, U = directions()
name = sys.argv[1]
tl = float(sys.argv[2]) if len(sys.argv) > 2 else 900
t0 = time.time()
if name.startswith("BD"):          # BDk: B_k(D) u conj(B_k(D))
    k = int(name[2:])
    B = sums_ball(D, k)
    P = unique_rows(np.concatenate([B, conj_rows(B)]))
    info = f"|B_{k}(D)| = {len(B)}, shared with conj: {2 * len(B) - len(P)}"
elif name.startswith("BU"):        # BUk: B_k(U)
    k = int(name[2:])
    P = sums_ball(U, k)
    info = f"|B_{k}(U)| = {len(P)}"
elif name.startswith("BDU"):
    raise SystemExit
elif name.startswith("B2D1U"):     # B_2(D) u conj(B_2(D)) u B_1(U)+B_1(U) ... not used
    raise SystemExit
else:
    raise SystemExit("unknown")
n = len(P)
E, J = build_edges(P, U)
S = PointSet(P)
i0 = S.index([[0] * 12])[0]; i1 = S.index([[7] + [0] * 11])[0]; i2 = S.index([D[14]])[0]
tri = (int(i0), int(i1), int(i2))
assert all(x >= 0 for x in tri)
Es = set(map(tuple, np.sort(E, axis=1).tolist()))
assert (min(i0, i1), max(i0, i1)) in Es and (min(i0, i2), max(i0, i2)) in Es and (min(i1, i2), max(i1, i2)) in Es
deg = np.bincount(E.ravel(), minlength=n)
inO = int(np.sum(np.all(P % 7 == 0, axis=1)))
print(f"{name}: {info}; vertices {n} (in Z[zeta21]: {inO}), edges {len(E)}, degree min/mean/max "
      f"{deg.min()}/{deg.mean():.1f}/{deg.max()}  [build {time.time() - t0:.1f}s]", flush=True)
cls = cnf_clauses(n, E.tolist(), tri)
s = Solver(name="cadical195", bootstrap_with=cls)
t1 = time.time()
r = solve_limited(s, tlimit=tl)
dt = time.time() - t1
if r is True:
    col = colouring_from_model(s.get_model(), n)
    print(f"{name}: SAT (4-colourable) in {dt:.1f}s; colouring verified: {check_colouring(col, E.tolist())}", flush=True)
    np.save(f"col_{name}.npy", np.array(col))
elif r is False:
    print(f"{name}: UNSAT in {dt:.1f}s (pysat, not certified)", flush=True)
else:
    print(f"{name}: UNKNOWN after {dt:.1f}s", flush=True)
np.save(f"pts_{name}.npy", P)
