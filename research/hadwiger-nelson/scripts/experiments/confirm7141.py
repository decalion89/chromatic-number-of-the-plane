"""Confirm the surviving pairs, one representative per orbit.

Forcing is equivariant: if (u,v) is forced then so is (g(u), g(v)) for every g
in the group, and if it is not then neither are they.  So the 2301 pairs that
survive the filter need only one solver call per ORBIT -- six times fewer --
and any forced pair found comes with five more for nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

d = json.load(open(HN_DIR + "/data/five_symmetric.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
g = build_graph(P); n = g.n; K = 5
pos = {p: i for i, p in enumerate(g.vertices)}
rot = _rot60(F)
perm = [pos[rot(g.vertices[v])] for v in range(n)]
cand = [tuple(x) for x in json.load(open(
    "/tmp/hn/cand7141.json"))]
seen, reps = set(), []
for a, b in cand:
    key = (min(a, b), max(a, b))
    if key in seen: continue
    x, y = a, b
    for _ in range(6):
        seen.add((min(x, y), max(x, y)))
        x, y = perm[x], perm[y]
    reps.append((a, b))
print(f"{len(cand)} candidates -> {len(reps)} orbit representatives", flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
t0 = time.time()
s = Solver(name="cd19", bootstrap_with=cnf)
forced = []
for i, (a, b) in enumerate(reps):
    if not s.solve(assumptions=[X(a, 0), X(b, 1)]):
        dd = (P[a] - P[b]).norm2()
        forced.append((a, b, str(dd)))
        print(f"  *** FORCED EQUAL: v{a} v{b}  d^2={dd}  ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
    if i < 5 or (i + 1) % 25 == 0:
        print(f"    {i+1}/{len(reps)}  forced so far {len(forced)}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
s.delete()
print(f"\n{len(forced)} forced orbits of {len(reps)}   [{time.time()-t0:.0f}s]",
      flush=True)
json.dump(forced, open("/tmp/hn/forced7141.json", "w"))
