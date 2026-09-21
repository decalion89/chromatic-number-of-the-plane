"""Wide against tight: measure the thing the chromatic number does not see.

Criticality arrives at 207 of Sa's 397 points and rigidity only at 397, so
something in that gap is doing the work and it is not the chromatic number.
Naming it is worth more than another search.

Two candidates are measurable on the same peel.  Saturation -- how many Moser
spindles a point belongs to -- is the account this repository has been giving,
and it should climb across the gap if it is the right one.  Mean degree is the
obvious rival and has already been shown decoupled from chi, so it should NOT
climb, or should climb far less.

Counting spindles exactly is expensive, but the spindle's defining feature is
cheap: a Moser spindle is built from two rhombi, and a rhombus is a pair of
points at squared distance 3 with two common unit neighbours.  So count, per
point, the pairs at squared distance 3 through it that close into rhombi --
the raw material the spindle is assembled from -- and watch it against size.
"""
import sys, time, math, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete

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
rng.shuffle(order)
print(f"Sa {len(P0)} points, same peel order as the census curve"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def measure(P):
    bb = IntBasis.covering(P)
    rr = bb.rows(P)
    assert bb.overflow_headroom(rr) < 1.0
    n = len(P)
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(bb, rr)))
    nb = [set() for _ in range(n)]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
    dm2, D22 = bb.dim, bb.D * bb.D
    rhombi = 0
    per = defaultdict(int)
    for i in range(n):
        dd = rr - rr[i]
        s = bb._field_square(dd[:, :dm2]) + bb._field_square(dd[:, dm2:])
        good = s[:, 0] == 3 * D22
        for m in range(1, dm2):
            good &= s[:, m] == 0
        for j in np.nonzero(good)[0]:
            j = int(j)
            if j <= i:
                continue
            shared = nb[i] & nb[j]
            if len(shared) >= 2:
                rhombi += 1
                per[i] += 1
                per[j] += 1
                for x in shared:
                    per[x] += 1
    return n, len(E), rhombi, (sum(per.values()) / n if n else 0)


for size in (7, 107, 157, 207, 257, 307, 347, 397):
    P = list(PIN) + order[:size - len(PIN)]
    n, ne, rh, sat = measure(P)
    print(f"   {n:4d} points, {ne:5d} edges, {ne/n:5.2f}/v: {rh:5d} rhombi, "
          f"{sat:6.2f} memberships per point  [{time.time()-t0:.0f}s]",
          flush=True)
print("DONE", flush=True)
