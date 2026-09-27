"""Use the group to multiply the colourings for free.

The filter needs many genuinely different colourings, and on the 7141-vertex
symmetric graph each one costs cadical minutes -- a full scan is out of reach,
and blocking a sample diversifies too weakly to make up for it: nine colourings
left five million of twenty-five million pairs alive.

But the graph is C6-invariant, and that is worth six colourings for the price
of one.  If c is proper and g is in the group then c . g is proper too, and it
agrees on (u,v) exactly when c agrees on (g(u), g(v)).  So a pair survives the
expanded family only if c identifies its WHOLE orbit -- a far stronger test --
and the expansion costs one permutation of the vertex indices, computed once.

Everything here is still a certificate in the only direction that matters: a
pair that differs in any of these colourings is proved not forced, because
each of them is a genuine proper 5-colouring of the graph.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

NCOL = int(sys.argv[1]) if len(sys.argv) > 1 else 6
d = json.load(open(HN_DIR + "/data/five_symmetric.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
g = build_graph(P); n = g.n; K = 5
pos = {p: i for i, p in enumerate(g.vertices)}
rot = _rot60(F)
perms = []
cur = list(range(n))
for _ in range(6):
    perms.append(cur)
    cur = [pos[rot(g.vertices[v])] for v in cur]
print(f"n={n}, group action as {len(perms)} vertex permutations "
      f"(identity check: {perms[0][:3]})", flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
for v in range(n):
    for c in range(K):
        for e in range(c + 1, K):
            cnf.append([-X(v, c), -X(v, e)])
t0 = time.time()
s = Solver(name="cd19", bootstrap_with=cnf)
if not s.solve():
    print("  *** IT REFUSES FIVE COLOURS ***", flush=True); sys.exit()
def readout():
    p = set(l for l in s.get_model() if l > 0)
    return [next(c for c in range(K) if X(v, c) in p) for v in range(n)]
raw = [readout()]
rng = random.Random(3)
def expand(cols):
    out = []
    for c in cols:
        for pm in perms:
            out.append([c[pm[v]] for v in range(n)])
    return out
def surviving(cols):
    b = defaultdict(int)
    for v in range(n): b[tuple(c[v] for c in cols)] += 1
    return sum(k * (k - 1) // 2 for k in b.values())
print(f"  1 raw -> {len(expand(raw))} effective: "
      f"{surviving(expand(raw))} surviving pairs of {n*(n-1)//2}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
while len(raw) < NCOL:
    s.add_clause([-X(v, raw[-1][v]) for v in rng.sample(range(n), 30)])
    if not s.solve(): break
    raw.append(readout())
    eff = expand(raw)
    print(f"  {len(raw)} raw -> {len(eff)} effective: {surviving(eff)} surviving"
          f"   [{time.time()-t0:.0f}s]", flush=True)
s.delete()
eff = expand(raw)
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in eff)].append(v)
cand = [(a, b) for vs in buck.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i+1:]]
print(f"\n  {len(raw)} raw colourings, {len(eff)} effective, "
      f"{len(cand)} candidate pairs   [{time.time()-t0:.0f}s]", flush=True)
json.dump([[int(a), int(b)] for a, b in cand[:500000]],
          open("/tmp/hn/cand7141.json", "w"))
