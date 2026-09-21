"""Reflect instead of rotate: the global operation nobody here has tried.

Every tightening operation used in this work is a rotation.  The bite unions a
graph with a rotated copy sharing the fixed POINT, and the matching between
them is the ring, whose points the rotation moves by exactly one.  A reflection
is the other kind of planar isometry, and it gives a different structure: the
two copies share the whole mirror LINE, and p meets its image at distance one
exactly when p lies at distance 1/2 from the mirror.

So the analogue of the ring is a PAIR OF LINES at distance 1/2 either side, and
the analogue of "how many points on the ring" is how many points sit on them.
That count is what to maximise, and it is cheap: for each direction the field
offers, project every point onto the normal and look for the offset with the
most points half a unit away on both sides.

Sa is already symmetric under reflection in the axes through its centre, so
those mirrors give nothing.  Every other line is new.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
t0 = time.time()
TEN = [
    [[0, 1, 2, 3, 4, 5, 6]],
    [[1, 2], [0, 3, 4, 5, 6]], [[2, 3], [0, 1, 4, 5, 6]],
    [[0, 1, 2, 3, 4], [5, 6]], [[3, 4], [0, 1, 2, 5, 6]],
    [[0, 1, 4], [2, 3, 5, 6]], [[0, 2, 3, 4, 5], [1, 6]],
    [[1, 2, 4, 5], [0, 3, 6]], [[4, 5], [0, 1, 2, 3, 6]],
    [[0, 2, 5], [1, 3, 4, 6]],
]
P0 = build_Sa(K)
import math
b0 = IntBasis.covering(P0)
r0 = b0.rows(P0)
dm, d2 = b0.dim, b0.D * b0.D
zi = P0.index(Point(K.zero(), K.zero()))
d = r0 - r0[zi]
sq = b0._field_square(d[:, :dm]) + b0._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring0 = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
ring0.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
SEVEN = [P0[zi]] + [P0[i] for i in ring0]
print(f"Sa: {len(P0)} points  [{time.time()-t0:.0f}s]", flush=True)

half = K.rational(Fr(1, 2))
# the directions the field gives as unit normals: multiples of thirty degrees
NORMALS = []
c30, s30 = K.sqrt(3) * half, half
rot30 = Rotation(c30, s30)
n = Point(K.rational(1), K.zero())
for _ in range(12):
    NORMALS.append(n)
    n = rot30(n)

# for each normal, the projections of every point, and the offsets that put
# the most points exactly half a unit either side
best, OFFSETS = [], []
for ni, nrm in enumerate(NORMALS):
    proj = [p.x * nrm.x + p.y * nrm.y for p in P0]
    cnt = Counter(proj)
    for t in set(proj):
        # the mirror at offset t + 1/4 pairs the plane t with the plane t+1/2
        lo, hi = t, t + half
        pair = min(cnt.get(lo, 0), cnt.get(hi, 0))
        if pair >= 4:
            # field elements are not orderable, so the offset is held by
            # index and never enters the sort key
            OFFSETS.append(t)
            best.append((pair, cnt.get(lo, 0) + cnt.get(hi, 0), ni,
                         len(OFFSETS) - 1))
best.sort(reverse=True)
print(f"{len(best)} (direction, offset) mirrors with 4+ matched points; "
      f"best {[b[0] for b in best[:8]]}  [{time.time()-t0:.0f}s]", flush=True)


def census(U, label):
    b = IntBasis.covering(U)
    r = b.rows(U)
    if b.overflow_headroom(r) >= 1.0:
        return None, 0
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    m = len(U)
    W = [U.index(q) for q in SEVEN]
    cls = [[1 + v * k + c for c in range(k)] for v in range(m)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        print(f"*** {label}: {m} points NOT 4-COLOURABLE ***", flush=True)
        return 0, len(E)
    keep = sum(1 for p in TEN
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(p)
                                        for e in blk]))
    sv.delete()
    return keep, len(E)


bestc, tried = 10, 0
for pair, tot, ni, ti in best[:60]:
    t = OFFSETS[ti]
    nrm = NORMALS[ni]
    mid = t + half * half              # the mirror plane, t + 1/4
    U, seen = list(P0), set(P0)
    for p in P0:
        s = p.x * nrm.x + p.y * nrm.y
        lam = (mid - s) * K.rational(2)
        q = Point(p.x + nrm.x * lam, p.y + nrm.y * lam)
        if q not in seen:
            seen.add(q)
            U.append(q)
    c, ne = census(U, f"mirror dir {ni} offset {t}")
    tried += 1
    if c is None:
        continue
    if c < bestc:
        bestc = c
        print(f"*** direction {ni}, offset {t}, {pair} matched: {len(U)} "
              f"points, {ne} edges, census {c} of 10 -- new best "
              f"(the bite gives 3) ***  [{time.time()-t0:.0f}s]", flush=True)
    if tried % 10 == 0:
        print(f"   {tried} mirrors, best census {bestc}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} mirrors, best census {bestc} of 10"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
