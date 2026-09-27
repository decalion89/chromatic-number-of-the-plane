"""Check directly that the pruning is what costs the other two pairs.

Y's census gives seven surviving patterns and four of them have the centre
alone, which the hand derivation says cannot happen in the unpruned union.
Since Y is the union LESS two vertices, and removing vertices can only add
survivors, the union's survivors are a subset of Y's -- so seven SAT calls
settle exactly which of the seven it keeps.

Predicted: the three with the centre beside an antipodal pair, and no others.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Sb, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
t0 = time.time()
SEVEN = [
    [[0, 1, 4], [2, 3, 5, 6]],
    [[1, 2, 4, 5], [0, 3, 6]],
    [[0], [1], [2, 4, 5], [3, 6]],
    [[0], [2], [1, 4, 5], [3, 6]],
    [[0, 2, 5], [1, 3, 4, 6]],
    [[0], [1, 2, 4], [5], [3, 6]],
    [[0], [4], [1, 2, 5], [3, 6]],
]
ZERO = Point(K.zero(), K.zero())


def seven_points(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    ci = P.index(ZERO)
    key = {tuple(r[i]): i for i in range(n)}
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
    ms = set(ring)
    pairs = [(j, key[tuple(2 * r[ci] - r[j])]) for j in ring
             if tuple(2 * r[ci] - r[j]) in key
             and key[tuple(2 * r[ci] - r[j])] in ms
             and key[tuple(2 * r[ci] - r[j])] > j]
    six = [x for pr in pairs[:3] for x in pr]
    six.sort(key=lambda i: math.atan2(float(P[i].y), float(P[i].x)))
    return n, E, [ci] + six


def survivors(P, label):
    n, E, W = seven_points(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    assert sv.solve(), f"{label} does not colour"
    keep = []
    for part in SEVEN:
        lits = [1 + W[e] * k + bi
                for bi, blk in enumerate(part) for e in blk]
        if sv.solve(assumptions=lits):
            keep.append(part)
    sv.delete()
    alone = sum(1 for p in keep if [0] in p)
    print(f"{label}: {n} pts, {len(E)} edges -> {len(keep)} of 7 survive, "
          f"centre alone in {alone}  [{time.time()-t0:.0f}s]", flush=True)
    for p in keep:
        print(f"    {p}", flush=True)
    return keep


survivors(build_Y(K), "Y (pruned)")
U, seen = [], set()
for p in build_Sa(K) + build_Sb(K):
    if p not in seen:
        seen.add(p)
        U.append(p)
survivors(U, "Sa u Sb (unpruned)")
print("DONE", flush=True)
