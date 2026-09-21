"""Enumerate the global operations and measure each one's tightening.

Rigidity is global: no local statistic predicts, produces or improves it, and
the only lever that has ever moved the census here is the bite -- a union with
a rotated copy sharing the centre the rotation fixes.  So enumerate the bites.

A bite is a pair (centre, ring): a point c of the graph and a rational squared
radius D about it whose joining rotation lies in the field, so that ring points
move by exactly one and the two copies meet along a matching.  De Grey used
c = origin, D = 4, and it takes Sa's census from ten to three.  Every other
pair is a different global operation and none has been measured.

Measured by TIGHTENING rather than by forcing, which is the weaker and more
sensitive question: adding the copy can only remove survivors, so testing Sa's
own ten patterns is exact and costs at most ten calls.  A census of zero would
mean the union does not 4-colour -- a 5-chromatic unit-distance graph in at
most 794 points, against G's 1581.
"""
import sys, time, math
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
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
print(f"Sa: {len(P0)} points; seven test points fixed  [{time.time()-t0:.0f}s]",
      flush=True)

# every (centre, ring) the field can bite, ordered by ring population
cands = []
for ci in range(len(P0)):
    dd = r0 - r0[ci]
    s = b0._field_square(dd[:, :dm]) + b0._field_square(dd[:, dm:])
    g = np.ones(len(s), dtype=bool)
    for m in range(1, dm):
        g &= s[:, m] == 0
    grp = defaultdict(int)
    for off in np.nonzero(g)[0]:
        v = Fr(int(s[off, 0]), d2)
        if v:
            grp[v] += 1
    for D, cnt in grp.items():
        if cnt >= 2 and closable_distance(D):
            cands.append((cnt, ci, D))
cands.sort(reverse=True)
seen_shape = set()
short = []
for cnt, ci, D in cands:
    key = (cnt, D, tuple(sorted(str(x) for x in (P0[ci].x.c[0],))))
    short.append((cnt, ci, D))
print(f"{len(short)} (centre, ring) bites available  [{time.time()-t0:.0f}s]",
      flush=True)

best = 10
tried = 0
for cnt, ci, D in short:
    rot = rotation_joining(D, K)
    for sgn in (+1, -1):
        f = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(P0[ci])
        s2 = set(P0)
        U = list(P0) + [q for q in (f(p) for p in P0) if q not in s2]
        b = IntBasis.covering(U)
        r = b.rows(U)
        if b.overflow_headroom(r) >= 1.0:
            continue
        E = sorted(set((min(a, c), max(a, c))
                       for a, c in fast_edges_complete(b, r)))
        n = len(U)
        W = [U.index(q) for q in SEVEN]
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, c in E:
            for col in range(k):
                cls.append([-(1 + a * k + col), -(1 + c * k + col)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        tried += 1
        if not sv.solve():
            print(f"*** centre {ci}, ring D={D}, dir {sgn}: {n} points NOT "
                  f"4-COLOURABLE -- 5-chromatic below G's 1581 ***"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            sv.delete()
            sys.exit()
        keep = sum(1 for p in TEN
                   if sv.solve(assumptions=[1 + W[e] * k + bi
                                            for bi, blk in enumerate(p)
                                            for e in blk]))
        sv.delete()
        if keep < best:
            best = keep
            print(f"*** centre {ci}, ring D={D} ({cnt} pts), dir {sgn}: "
                  f"{n} points, census {keep} of 10 -- new best ***"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
        if tried % 10 == 0:
            print(f"   {tried} bites, best census {best}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} bites, best census {best} of 10  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
