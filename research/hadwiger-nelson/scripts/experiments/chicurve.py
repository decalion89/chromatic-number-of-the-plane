"""Does the census fall exactly when the graph reaches its chromatic number?

The census sits at the ceiling to 307 of Sa's 397 points and then falls off a
cliff.  The natural explanation is that four colours are simply not scarce
below that size: if a 307-point subgraph is 3-COLOURABLE then its 4-colourings
are abundant and of course every pattern survives, and rigidity would appear
exactly when the fourth colour becomes necessary.

That is a claim with a cheap test.  Walk the same peeling order and ask for
the chromatic number at each size -- 3-colourable or not -- alongside the
census reading already measured.  If the two thresholds coincide, the census
is a proxy for criticality and nothing more; if the graph is 4-chromatic well
before the census moves, they are different phenomena and the gap between them
is the thing to understand.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
rng = random.Random(11)
P0 = build_Sa(K)
b = IntBasis.covering(P0)
r = b.rows(P0)
dm, d2 = b.dim, b.D * b.D
zi = P0.index(Point(K.zero(), K.zero()))
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
ring.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
PIN = [P0[zi]] + [P0[i] for i in ring]
order = [p for p in P0 if p not in set(PIN)]
rng.shuffle(order)                      # the SAME order as the census curve
print(f"Sa {len(P0)} points; seven pinned, same peel order as the census"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def chi_at_most(P, k):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    out = s.solve()
    s.delete()
    return out, len(E)


CENSUS4 = {7: 715, 27: 715, 57: 715, 107: 715, 157: 715, 207: 715,
           257: 715, 307: 715, 347: 577, 397: 10}
for size in sorted(CENSUS4):
    P = list(PIN) + order[:size - len(PIN)]
    three, ne = chi_at_most(P, 3)
    four, _ = chi_at_most(P, 4)
    chi = 3 if three else (4 if four else 5)
    print(f"   {len(P):4d} points, {ne:5d} edges: chi = {chi}, "
          f"census at four = {CENSUS4[size]} of 715"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
