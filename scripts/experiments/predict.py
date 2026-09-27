"""Use the saturation criterion to predict, instead of only to explain.

Across the gap between criticality and rigidity, density rises by a factor of
two and rhombus saturation by a factor of eight, which says saturation is what
matters.  That is an explanation, and an explanation that cannot predict is
worth little.  So predict: at a size where a RANDOM subset of Sa reads the
ceiling -- 307 points, saturation 1.77, census 715 of 715 -- a subset of the
same size CHOSEN for saturation should read lower, and if the criterion is any
good it should read much lower.

Chosen greedily: keep the seven pinned points, then repeatedly add whichever
remaining point of Sa joins the most rhombi with what is already there.  Sizes
are matched to the random curve exactly so the two columns compare.

If a 307-point subset censuses near ten it is also a smaller rigid graph than
Sa, which would be worth having on its own.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
P0 = build_Sa(K)
n0 = len(P0)
b = IntBasis.covering(P0)
r = b.rows(P0)
dm, d2 = b.dim, b.D * b.D
E0 = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
nb = [set() for _ in range(n0)]
for a, c in E0:
    nb[a].add(c)
    nb[c].add(a)
zi = P0.index(Point(K.zero(), K.zero()))
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
ring.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
PIN = [zi] + ring

# every rhombus of Sa, once
RH = []
for i in range(n0):
    dv = r - r[i]
    s = b._field_square(dv[:, :dm]) + b._field_square(dv[:, dm:])
    good = s[:, 0] == 3 * b.D * b.D
    for m in range(1, dm):
        good &= s[:, m] == 0
    for j in np.nonzero(good)[0]:
        j = int(j)
        if j <= i:
            continue
        shared = sorted(nb[i] & nb[j])
        if len(shared) >= 2:
            for a in range(len(shared) - 1):
                for c in range(a + 1, len(shared)):
                    RH.append((i, j, shared[a], shared[c]))
print(f"Sa: {n0} points, {len(E0)} edges, {len(RH)} rhombi (as 4-sets)"
      f"  [{time.time()-t0:.0f}s]", flush=True)

member = defaultdict(list)
for idx, q in enumerate(RH):
    for v in q:
        member[v].append(idx)

chosen = set(PIN)
complete = set(idx for idx, q in enumerate(RH) if set(q) <= chosen)
covered = defaultdict(int)
order = []
while len(chosen) < n0:
    best, bestgain = None, -1
    for v in range(n0):
        if v in chosen:
            continue
        gain = sum(1 for idx in member[v]
                   if idx not in complete and
                   len(set(RH[idx]) - chosen - {v}) == 0)
        if gain > bestgain:
            best, bestgain = v, gain
    chosen.add(best)
    order.append(best)
    for idx in member[best]:
        if idx not in complete and set(RH[idx]) <= chosen:
            complete.add(idx)


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
k = 4
PARTS = [p for p in ALL if len(p) <= k]
print(f"greedy order built  [{time.time()-t0:.0f}s]", flush=True)

for size in (207, 257, 307, 347, 397):
    idxs = PIN + order[:size - len(PIN)]
    P = [P0[i] for i in idxs]
    bb = IntBasis.covering(P)
    rr = bb.rows(P)
    assert bb.overflow_headroom(rr) < 1.0
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(bb, rr)))
    S = set(idxs)
    rh = sum(1 for q in RH if set(q) <= S)
    per = defaultdict(int)
    for q in RH:
        if set(q) <= S:
            for v in q:
                per[v] += 1
    sat = sum(per.values()) / len(P)
    W = list(range(7))                      # PIN comes first in idxs
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"   {len(P):4d} points: NOT 4-COLOURABLE", flush=True)
        sv.delete()
        continue
    surv = sum(1 for part in PARTS
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(part)
                                        for e in blk]))
    sv.delete()
    print(f"   {len(P):4d} points, {len(E):5d} edges ({len(E)/len(P):.2f}/v), "
          f"{rh:5d} rhombi, {sat:5.2f} per point: census {surv} of 715"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
