"""Which unit does the rigidity threshold travel in -- gadgets, or incidences?

The overlap invariant turns the whole question into one number, and the number
depends on a choice that cannot be settled by argument.  At four colours the
census leaves the ceiling at 0.161 Moser spindles per point, which is 1.13
vertex incidences per point since a spindle has seven vertices.  Carried to
five colours with a 509-vertex gadget the two readings diverge by 143:

    threshold in GADGETS      0.161 per point  ->  6 new points per copy
                                                   98.8 per cent overlap
    threshold in INCIDENCES   1.13  per point  ->  127 new points per copy
                                                   75 per cent overlap, which
                                                   is Sa's own

The first says the construction is hopeless, the second that reproducing Sa's
own overlap proportion would do.  A third data point decides which travels.

Three colours, where the critical gadget is the TRIANGLE -- three vertices,
3-chromatic, and the triangular lattice is saturated with them.  Peel a lattice
patch and find where its 3-colourings stop being rigid, in both units.  If the
gadget count matches four colours' 0.161 the count travels; if the incidence
count matches 1.13 the incidences do.
"""
import sys, time, random, math
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
rng = random.Random(11)
half = K.rational(Fr(1, 2))
s3 = K.sqrt(3) * half
pts, seen = [], set()
R = 7
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        x = K.rational(a) + K.rational(b) * half
        y = K.rational(b) * s3
        p = Point(x, y)
        if float(x) ** 2 + float(y) ** 2 <= R * R and p not in seen:
            seen.add(p)
            pts.append(p)
print(f"triangular lattice patch: {len(pts)} points  [{time.time()-t0:.0f}s]",
      flush=True)


def stats(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    n = len(P)
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    nb = [set() for _ in range(n)]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
    tri = sum(1 for a, c in E for x in nb[a] & nb[c] if x > c)
    return n, E, tri


# seven points to census: a centre and six at squared distance 4 about it,
# which in this lattice is the hexagon of radius two
n0, E0, tri0 = stats(pts)
b = IntBasis.covering(pts)
r = b.rows(pts)
dm, d2 = b.dim, b.D * b.D
zi = pts.index(Point(K.zero(), K.zero()))
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
ring.sort(key=lambda i: math.atan2(float(pts[i].y), float(pts[i].x)))
assert len(ring) >= 6, len(ring)
ring = ring[:6]
PIN = [pts[zi]] + [pts[i] for i in ring]
order = [p for p in pts if p not in set(PIN)]
rng.shuffle(order)
print(f"{n0} points, {len(E0)} edges, {tri0} triangles "
      f"({tri0/n0:.3f}/pt, {3*tri0/n0:.3f} incidences/pt); ring of "
      f"{len(ring)}  [{time.time()-t0:.0f}s]", flush=True)


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


k = 3
PARTS = [p for p in partitions(list(range(7))) if len(p) <= k]
print(f"ceiling at three colours: {len(PARTS)} patterns", flush=True)
sizes = sorted({7, 20, 40, 60, 80, 100, 130, len(pts)})
for size in sizes:
    if size > len(pts):
        continue
    P = list(PIN) + order[:size - len(PIN)]
    n, E, tri = stats(P)
    W = list(range(7))
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"   {n:4d} points: NOT 3-colourable", flush=True)
        sv.delete()
        continue
    surv = sum(1 for part in PARTS
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(part)
                                        for e in blk]))
    sv.delete()
    mark = "  <- off the ceiling" if surv < len(PARTS) else ""
    print(f"   {n:4d} points, {len(E):5d} edges, {tri:5d} triangles "
          f"({tri/n:6.3f}/pt, {3*tri/n:6.3f} inc/pt): census {surv} of "
          f"{len(PARTS)}{mark}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
