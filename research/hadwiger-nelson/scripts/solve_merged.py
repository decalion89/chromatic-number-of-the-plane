"""Settle a saved gadget candidate: is H with a and b merged 5-colourable?

Merging a and b is the same question as c(a) = c(b), asked as a graph.  A
triangle of the merged graph is pinned to colours 0, 1, 2 (sound: any colouring
can be permuted into that), which cuts the search by up to 60.  Runs one solver
named on the command line with an optional conflict budget; SAT prints the
verdict and saves nothing, UNSAT points at scripts/verify_gadget.py.
"""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
t0 = time.time(); K = 5
path, name = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "cd19")
budget = int(sys.argv[3]) if len(sys.argv) > 3 else 0
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
a, b = d["A"], d["B"]
g = build_graph(P); n = g.n
ren = lambda v: a if v == b else (v if v < b else v - 1)
ME = sorted({(min(ren(x), ren(y)), max(ren(x), ren(y))) for x, y in g.edges()})
m = n - 1
nb = {}
for x, y in ME: nb.setdefault(x, set()).add(y); nb.setdefault(y, set()).add(x)
tri = next((x, y, z) for x, y in ME for z in nb[x] & nb[y])
X = lambda v, c: 1 + v * K + c
s = Solver(name=name)
for v in range(m): s.add_clause([X(v, c) for c in range(K)])
for x, y in ME:
    for c in range(K): s.add_clause([-X(x, c), -X(y, c)])
for k, v in enumerate(tri): s.add_clause([X(v, k)])
print(f"{path}: n={n}, merged graph {m} vertices {len(ME)} edges; solver {name}, budget {budget or 'none'}", flush=True)
if budget:
    s.conf_budget(budget); r = s.solve_limited()
else:
    r = s.solve()
print(f"  merged graph 5-colourable: {r}   conflicts {s.accum_stats().get('conflicts')}   [{time.time()-t0:.0f}s]", flush=True)
if r is False:
    print("  UNSAT -- run scripts/verify_gadget.py on this file before believing anything")
