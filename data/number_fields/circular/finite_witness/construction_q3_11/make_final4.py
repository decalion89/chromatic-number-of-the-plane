"""make_final4.py W.json NAME: the witness file NAME.json.gz (denominator, points, all unit pairs as edges, a proper
4-colouring, the cycles, the fixed vertex = the origin) and NAME.cnf (the formula of check_witness4.py)."""
import sys, json, gzip, itertools
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from check_witness4 import cnf_text, unit
from pysat.solvers import Solver
import numpy as np
W = json.load(open(sys.argv[1])); name = sys.argv[2]
D = W["D"]; P = [list(map(int, p)) for p in W["points"]]; n = len(P)
A = np.array(P, dtype=np.int64)
r = (A[:, None, :4] - A[None, :, :4]); s = (A[:, None, 4:] - A[None, :, 4:])
w = np.array([1, 3, 11, 33], dtype=np.int64)
r0 = (r * r * w).sum(-1) + (s * s * w).sum(-1)
E = []
for i, j in zip(*np.nonzero(r0 == D * D)):
    if i < j and unit(P[i], P[j], D):
        E.append([int(i), int(j)])
fixed = P.index([0] * 8)
cycles = [c for c in W["cycles"] if len(c) % 4 == 0]
x = lambda v, k: 4 * v + k + 1
sol = Solver(name="cadical153")
for v in range(n):
    sol.add_clause([x(v, k) for k in range(4)])
    for k, l in itertools.combinations(range(4), 2):
        sol.add_clause([-x(v, k), -x(v, l)])
for i, j in E:
    for k in range(4):
        sol.add_clause([-x(i, k), -x(j, k)])
assert sol.solve()
m = set(l for l in sol.get_model() if l > 0)
col = [next(k for k in range(4) if x(v, k) in m) for v in range(n)]
out = {"field": "Q(sqrt3, sqrt11)", "basis": "(1, sqrt3, sqrt11, sqrt33)", "denominator": D, "points": P,
       "edges": E, "colouring": col, "cycles": cycles, "fixed_vertex": fixed}
with gzip.open(name + ".json.gz", "wt") as f:
    json.dump(out, f)
open(name + ".cnf", "w").write(cnf_text(n, E, cycles, fixed))
print(f"{name}: {n} points, {len(E)} edges, {len(cycles)} cycles, fixed vertex {fixed}")
