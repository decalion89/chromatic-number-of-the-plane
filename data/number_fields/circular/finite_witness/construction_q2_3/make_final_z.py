"""make_final_z.py W.json NAME: the witness file NAME.json.gz over Q(sqrt2, sqrt3) and NAME.cnf (the formula of
check_witness4.py). W.json holds points of Q(zeta_24) (coordinates over the power basis 1, zeta, ..., zeta^7, over
the denominator D); a point z = x + iy becomes the plane point (x, y) with coordinates over (1, sqrt2, sqrt3, sqrt6)
and denominator 4D (to_xy.py). Edges: all unit pairs, exactly. A proper 4-colouring from CaDiCaL. The fixed vertex
is the origin."""
import sys, os, json, gzip, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                  # to_xy.py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # check_witness4.py
from to_xy import xy
from check_witness4 import cnf_text, unit
from pysat.solvers import Solver
import numpy as np
W = json.load(open(sys.argv[1])); name = sys.argv[2]
D = 4 * W["D"]
P = [xy(c) for c in W["points"]]; n = len(P)
if len(set(map(tuple, P))) != n:
    sys.exit("repeated point")
A = np.array(P, dtype=np.int64)
w = np.array([1, 2, 3, 6], dtype=np.int64)
E = []
for i in range(n):                                   # rational part of the squared length must be D^2
    d = A[i + 1:] - A[i]
    r0 = (d[:, :4] * d[:, :4] * w).sum(1) + (d[:, 4:] * d[:, 4:] * w).sum(1)
    for j in np.nonzero(r0 == D * D)[0]:
        j = i + 1 + int(j)
        if unit(P[i], P[j], D, 2, 3):
            E.append([i, j])
fixed = P.index([0] * 8)
cycles = [c for c in W["cycles"] if len(c) % 4 == 0]
Es = set(map(tuple, E))
for c in cycles:
    for t in range(len(c)):
        a, b = c[t], c[(t + 1) % len(c)]
        if (min(a, b), max(a, b)) not in Es:
            sys.exit(f"cycle {c} uses a non-edge")
x = lambda v, k: 4 * v + k + 1
sol = Solver(name="cadical153")
for v in range(n):
    sol.add_clause([x(v, k) for k in range(4)])
    for k, l in itertools.combinations(range(4), 2):
        sol.add_clause([-x(v, k), -x(v, l)])
for i, j in E:
    for k in range(4):
        sol.add_clause([-x(i, k), -x(j, k)])
if not sol.solve():
    sys.exit("not 4-colourable")
m = set(l for l in sol.get_model() if l > 0)
col = [next(k for k in range(4) if x(v, k) in m) for v in range(n)]
out = {"field": "Q(sqrt2, sqrt3)", "basis": "(1, sqrt2, sqrt3, sqrt6)", "denominator": D, "points": P,
       "edges": E, "colouring": col, "cycles": cycles, "fixed_vertex": fixed}
with open(name + ".json.gz", "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as f:   # reproducible bytes
    f.write(json.dumps(out).encode())
open(name + ".cnf", "w").write(cnf_text(n, E, cycles, fixed))
print(f"{name}: {n} points, {len(E)} edges, {len(cycles)} cycles, fixed vertex {fixed}")
