"""How fast does the census fall with size, and where would five need to fall?

At four colours Sa scores 10 of 715 and every graph at five colours scores 855
of 855, up to the 13356-point union of ten translates.  Two endpoints are not
a curve.  Peeling Sa down and censusing as it shrinks gives the whole shape at
four, and the size at which it leaves the ceiling is the quantity to compare
against: at five nothing has left the ceiling at any size yet measured, so the
ratio of the two sizes bounds how much further five would have to go.

Vertices are dropped in random order with the seven pinned, so the curve is a
property of size rather than of any clever choice, and the same order is used
at both k so the two curves are comparable point by point.
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
ZERO = Point(K.zero(), K.zero())
b = IntBasis.covering(P0)
r = b.rows(P0)
dm, d2 = b.dim, b.D * b.D
zi = P0.index(ZERO)
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
ring.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
PIN = [P0[zi]] + [P0[i] for i in ring]
print(f"Sa {len(P0)} points, seven pinned  [{time.time()-t0:.0f}s]",
      flush=True)


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


ALL = list(partitions(list(range(7))))
order = [p for p in P0 if p not in set(PIN)]
rng.shuffle(order)


def census(P, k):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    W = [P.index(q) for q in PIN]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None, len(E)
    parts = [p for p in ALL if len(p) <= k]
    surv = sum(1 for part in parts
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(part)
                                        for e in blk]))
    sv.delete()
    return surv, len(E)


sizes = [len(PIN) + m for m in (0, 20, 50, 100, 150, 200, 250, 300, 340,
                                len(order))]
for k in (4, 5):
    tot = sum(1 for p in ALL if len(p) <= k)
    print(f"\nk={k}, ceiling {tot}", flush=True)
    for s in sizes:
        P = list(PIN) + order[:s - len(PIN)]
        surv, ne = census(P, k)
        tag = "UNCOLOURABLE" if surv is None else (
            f"{surv} of {tot}" + ("  <- off the ceiling" if surv < tot else ""))
        print(f"   {len(P):4d} points, {ne:5d} edges: {tag}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
