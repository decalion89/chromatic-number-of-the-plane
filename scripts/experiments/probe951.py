"""The smaller 5-chromatic graph, measured the same way as the others.

Overlap is what defeats the sigma-argument, so the question for any 5-chromatic
carrier is how much of itself a glue can reuse.  On Sa at four colours the best
glue reuses 56% and forcing appears; on Z at five the best rot60 about a vertex
reused 40% and nothing was forced.  This asks the same of the 951-vertex graph,
and measures its slack while it is at it.
"""
import os, sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, rotation_joining
from hn.graph import build_graph
from pysat.solvers import Solver

d = json.load(open("data/five_247_b.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
K = 5
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf)
s.solve()
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
r = random.Random(4)
for _ in range(60):
    if len(cols) >= 24: break
    s.set_phases([(1 if r.random() < 0.2 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): continue
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
free = min(sum(1 for v in range(n) if len({col[u] for u in g.adj[v]} | {col[v]}) < K)
           for col in cols)
buck = defaultdict(list)
for v in range(n):
    buck[tuple(c[v] for c in cols)].append(v)
cand = sum(len(b) * (len(b) - 1) // 2 for b in buck.values())
print(f"951-vertex graph: n={n} m={m} deg={2.0*m/n:.2f}  free@5={100.0*free/n:.2f}%"
      f"  forced-pair candidates at 5: {cand}", flush=True)

S = set(P)
t0 = time.time()
best = []
for name, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3)), ("rot180", Fr(1, 4))):
    r0 = rotation_joining(d2, F)
    for c in P:
        rot = r0.about(c)
        ov = sum(1 for p in P if rot(p) in S)
        if 1 < ov < n: best.append((ov, name))
best.sort(reverse=True)
print(f"  best two-copy glue overlaps: {[b[0] for b in best[:8]]} of {n}"
      f"  ({100.0*best[0][0]/n:.0f}% at most)   [{time.time()-t0:.0f}s]", flush=True)
