"""Does any graph here FORCE a distance at five colours?  A disjunction, not a pair.

Every forcing test in this project has asked the strongest possible question --
is THIS pair monochromatic in every 5-colouring -- and the pricing showed why
the answer is always no: a forced pair needs a graph already one vertex short of
6-chromatic.  The weak question has never been asked:

    is there a distance d such that EVERY 5-colouring of H has SOME
    monochromatic pair at distance d ?

That is a disjunction over all the d-pairs at once, so it is enormously weaker
than any single forced pair, and it is exactly the shape de Grey's argument
needed at four colours: a statement that some member of a finite list must
happen, consumed afterwards by copies.

It is also one solve.  "No monochromatic pair at distance d" is just the
d-pairs added as edges, so the question is whether H is still 5-colourable once
BOTH distances are forbidden:

    chi(H; {1, d}) > 5  <=>  H forces d at five colours.

Cheap, and never tried.  Asked of every squared distance H realises, richest
first, since a distance with more pairs constrains more.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
TOP = int(sys.argv[2]) if len(sys.argv) > 2 else 60
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
print(f"{NAME} n={n} unit edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

# every realised squared distance, bucketed by a rounded float and confirmed by
# exact arithmetic inside each bucket
buckets = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if dd > 16.0: continue
        buckets[round(dd, 7)].append((i, j))
print(f"  {len(buckets)} distance classes within d^2 <= 16   "
      f"[{time.time()-t0:.0f}s]", flush=True)
classes = sorted(buckets.items(), key=lambda kv: -len(kv[1]))
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])

hits = []
for rank, (d2, pairs) in enumerate(classes[:TOP]):
    if abs(d2 - 1.0) < 1e-9: continue
    exact = (g.vertices[pairs[0][0]] - g.vertices[pairs[0][1]]).norm2()
    pairs = [(i, j) for i, j in pairs
             if (g.vertices[i] - g.vertices[j]).norm2() == exact]
    cnf = list(base)
    for i, j in pairs:
        for c in range(K):
            cnf.append([-X(i, c), -X(j, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    ok = s.solve(); s.delete()
    tag = "" if ok else "   *** FORCES THIS DISTANCE AT FIVE ***"
    print(f"  d^2={d2:<12.6f} pairs={len(pairs):<6d} chi<=5 with both "
          f"forbidden: {ok}{tag}   [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        hits.append({"d2": d2, "pairs": len(pairs)})
json.dump({"graph": NAME, "forced_distances": hits},
          open(f"{ROOT}/data/two_distance_{NAME}", "w"))
print(f"\n  {len(hits)} forced distances   [{time.time()-t0:.0f}s]", flush=True)
