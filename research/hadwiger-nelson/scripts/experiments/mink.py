"""Rotations do not interlock.  Translations by unit vectors do.

The scan over C3, C6, C12 and the spindle group about six centres put free@5
between 14.2 % and 16.4 % -- every one LOOSER than the 803-point carrier's own
11.5 %.  Rotating a sprawling object about a point just lays copies beside each
other; the mean degree went 10.12 -> 10.79 for four times the vertices.

The arithmetic says why that can never work.  A vertex is free at k when some
colour is missing from its neighbourhood, and for a random-looking graph

    free@k  ~  (k-1) * (1 - 1/(k-1))^deg

At four colours Sa has degree 9.94, so the law predicts 3*(2/3)^9.94 = 5.2 %,
and Sa measures 0.00 % -- twenty times better than the law.  At five colours
the 803-point graph has degree 10.12, the law predicts 4*(3/4)^10.12 = 21.8 %,
and it measures 11.5 % -- under two times better.  Sa is STRUCTURALLY rigid;
our 5-chromatic graphs are merely random graphs that happen to refuse four.

The reason is where the chromatic number lives.  Sa refuses three everywhere at
once: every vertex sits deep inside several overlapping copies of the 39-point
diamond.  Z refuses four because of ONE spindle; strip the pivot and the rest
is 4-chromatic material, which at five colours is slack by construction.

So the 5-chromaticity has to be made local, and the way to do that is to cover
the object with overlapping copies of itself -- by TRANSLATION.  Translating by
a unit vector is free density: p and p+t are at distance exactly 1, so every
vertex gains an edge for every translation, before any coincidence.  H is the
hexagon plus its centre, the seven points of the standard unit-distance
7-set, and Z (+) H is the Minkowski sum: 7 copies of Z, each pair of them
offset by a hexagon edge, so the copies are threaded together rather than laid
side by side.

The hexagonal unit vectors need only sqrt3, which the field already has, so
this costs no new radical -- the graph stays in Q(sqrt3, sqrt11, sqrt247).
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
F = Field((3, 11, 247))
one, zero, half = F.rational(Fr(1)), F.rational(Fr(0)), F.rational(Fr(1, 2))
rt3_2 = F.sqrt(3) * half
HEX = [Point(zero, zero),
       Point(one, zero), Point(half, rt3_2), Point(-half, rt3_2),
       Point(-one, zero), Point(-half, -rt3_2), Point(half, -rt3_2)]

d = json.load(open(ROOT + "/data/five_247_c.json"))
Z = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]

def law(deg, K):
    return (K - 1) * ((K - 2) / (K - 1)) ** deg

def free_at(g, col, K):
    fr = 0
    for v in range(g.n):
        if len(set(col[u] for u in g.adj[v])) < K - 1:
            fr += 1
    return fr / g.n

def colour(g, K, budget=20_000_000):
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cl)
    s.conf_budget(budget)
    r = s.solve_limited()
    if r is not True:
        s.delete(); return None, r
    pos = set(l for l in s.get_model() if l > 0); s.delete()
    return [next(c for c in range(K) if X(v, c) in pos) for v in range(n)], True

def report(tag, V, K=5):
    g = build_graph(V); n = g.n; m = sum(len(a) for a in g.adj) // 2
    deg = 2.0 * m / n
    col, st = colour(g, K)
    if col is None:
        print(f"  {tag:22s} n={n} m={m} deg={deg:.2f}   *** {K}-colour "
              f"{'UNSAT' if st is False else 'budget out'} ***"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        if st is False:
            json.dump({"field_generators": list(F.gens), "n": n, "m": m,
                       "how": tag,
                       "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                                   [[c.numerator, c.denominator] for c in p.y.c]]
                                  for p in g.vertices]},
                      open(ROOT + "/data/six_candidate.json", "w"))
            print("  *** written data/six_candidate.json ***", flush=True)
        return g, None
    f = free_at(g, col, K)
    print(f"  {tag:22s} n={n} m={m} deg={deg:.2f}  free@{K}={100*f:6.3f}%"
          f"  (law {100*law(deg, K):6.3f}%, {law(deg, K)/max(f,1e-9):5.1f}x better)"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    return g, f

def sumset(A, B):
    seen, out = set(), []
    for t in B:
        for p in A:
            q = Point(p.x + t.x, p.y + t.y)
            if q not in seen:
                seen.add(q); out.append(q)
    return out

report("Z", Z)
Z1 = sumset(Z, HEX)
print(f"  Z (+) H : {len(Z1)} points   [{time.time()-t0:.0f}s]", flush=True)
report("Z (+) H", Z1)

# One more layer, if it is affordable.  Three copies of the hexagon rather than
# the full seven keeps the growth down while still threading the object.
TRI = [HEX[0], HEX[1], HEX[3], HEX[5]]
Z2 = sumset(Z1, TRI)
print(f"  Z (+) H (+) T : {len(Z2)} points   [{time.time()-t0:.0f}s]", flush=True)
if len(Z2) <= 26000:
    g2, f2 = report("Z (+) H (+) T", Z2)
    if f2 is not None and f2 < 0.02:
        json.dump({"field_generators": list(F.gens), "n": g2.n,
                   "m": sum(len(a) for a in g2.adj) // 2, "free_at_5": f2,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in g2.vertices]},
                  open(ROOT + "/data/tight_five.json", "w"))
        print("  written data/tight_five.json", flush=True)
