"""Does the symmetric 5-chromatic graph have a forced-equal pair at five?

This is the whole question.  A 5-chromatic unit-distance graph with two
vertices that every proper 5-colouring must give the same colour, at a distance
d with 4d^2 - 1 a square, is one rotation away from chi >= 6 -- that is exactly
de Grey's own step, verified on Y at four colours.

The graph here is the C6-invariant 7141-vertex spindle, the hardest instance
the project has built: cadical needs 322 seconds for one 5-colouring where a
6607-vertex graph of the same degree needs 4.  So the filter is expensive --
every colouring costs minutes rather than milliseconds -- but it is complete in
the only direction that matters: a pair that differs in any sampled colouring
is PROVED not forced, so the negatives are certificates and not samples.

Diversity comes from blocking a random sample of the last solution, which is
what works with cadical; phases would be ignored by it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

NCOL = int(sys.argv[1]) if len(sys.argv) > 1 else 12
d = json.load(open(HN_DIR + "/data/five_symmetric.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
g = build_graph(P); n = g.n; K = 5
print(f"n={n} m={sum(len(a) for a in g.adj)//2}, k={K}", flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
for v in range(n):                      # at-most-one: the colours are read off
    for c in range(K):
        for e in range(c + 1, K):
            cnf.append([-X(v, c), -X(v, e)])
t0 = time.time()
s = Solver(name="cd19", bootstrap_with=cnf)
if not s.solve():
    print("  *** IT REFUSES FIVE COLOURS ***", flush=True); sys.exit()
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
print(f"  colouring 1   [{time.time()-t0:.0f}s]", flush=True)
rng = random.Random(3)
while len(cols) < NCOL:
    last = cols[-1]
    s.add_clause([-X(v, last[v]) for v in rng.sample(range(n), 200)])
    if not s.solve(): break
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
    buck = defaultdict(int)
    for v in range(n): buck[tuple(c[v] for c in cols)] += 1
    cand = sum(k * (k - 1) // 2 for k in buck.values())
    print(f"  colouring {len(cols)}: {cand} surviving pairs of {n*(n-1)//2}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
s.delete()
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in cols)].append(v)
cand = [(a, b) for vs in buck.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i+1:]]
print(f"\n  {len(cols)} colourings, {len(cand)} candidates   [{time.time()-t0:.0f}s]",
      flush=True)
json.dump([[int(a), int(b)] for a, b in cand[:200000]],
          open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/cand7141.json", "w"))
