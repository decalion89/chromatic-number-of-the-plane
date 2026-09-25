"""What is the origin's relation to the ring it caps?

Sa's D = 4 ring takes at most two colours in every 4-colouring.  That is the
lemma, but it says nothing about WHICH two, and the obvious candidate for one
of them is the centre's own colour -- the ring sits at radius two from a
degree-thirty hub, close enough to be constrained by it and too far to be
adjacent.

The seven points are enumerable.  Force the centre and its six ring points to
each partition of seven things into at most k blocks and ask which survive.
What comes back is the exact joint statement, from which the cap is a
projection, and it says directly whether the centre is inside the ring's
classes or outside them.

Run at four colours, where the cap holds, and at five, where it does not, so
the difference is visible rather than asserted.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
P = build_Sa(K)
n = len(P)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
dm, d2 = b.dim, b.D * b.D
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
zi = P.index(Point(K.zero(), K.zero()))
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
import math
ring.sort(key=lambda i: math.atan2(float(P[i].y), float(P[i].x)))
W = [zi] + ring            # position 0 is the centre
print(f"Sa: {n} pts, {len(E)} edges; centre {zi}, ring {ring}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


for k in (4, 5):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    surv, total = [], 0
    for part in partitions(list(range(7))):
        if len(part) > k:
            continue
        total += 1
        lits = []
        for bi, blk in enumerate(part):
            for e in blk:
                lits.append(1 + W[e] * k + bi)
        if sv.solve(assumptions=lits):
            surv.append(tuple(tuple(sorted(x)) for x in part))
    sv.delete()
    withcentre = [p for p in surv
                  if any(0 in blk and len(blk) > 1 for blk in p)]
    alone = len(surv) - len(withcentre)
    ringcols = defaultdict(int)
    for p in surv:
        used = {bi for bi, blk in enumerate(p) for e in blk if e > 0}
        ringcols[len(used)] += 1
    print(f"\nk={k}: {len(surv)} of {total} patterns on centre+ring survive"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"   centre shares its colour with a ring point in {len(withcentre)};"
          f" centre alone in {alone}", flush=True)
    print(f"   colours used ON THE RING: {dict(sorted(ringcols.items()))}",
          flush=True)
    for p in surv[:10]:
        tag = "centre with ring" if any(0 in blk and len(blk) > 1
                                        for blk in p) else "centre alone"
        print(f"      {p}   [{tag}]", flush=True)
    if len(surv) > 10:
        print(f"       ... and {len(surv)-10} more", flush=True)
print("DONE", flush=True)
