"""What does Sa actually say about its hexagon?

The bite adds almost nothing.  Y is Sa together with one rotated copy, and the
two share ONE vertex and exactly SIX edges -- the matching that sends each
hexagon point to its image, which the rotation puts at distance one.  Six
edges, and the union forces a named pair at four colours where neither half
forces anything.  "One of three antipodal pairs is monochromatic" cannot do
that: it is consistent with the two halves naming different pairs.

So Sa says more than that about its hexagon, and the exact statement is
enumerable.  The six points of the D = 4 ring carry no edges among themselves
-- consecutive ones are two apart, antipodal ones four -- so a priori every
partition of six things into at most four blocks is available.  Ask the solver
which ones actually occur: force the hexagon to each partition in turn and see
whether Sa still colours.  What survives IS the lemma, stated exactly, with no
paraphrase in the way.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from itertools import product
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
Sa = build_Sa(K)
ZERO = Point(K.zero(), K.zero())
b = IntBasis.covering(Sa)
r = b.rows(Sa)
dm, d2 = b.dim, b.D * b.D
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(Sa)
zi = Sa.index(ZERO)
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
# order the hexagon by angle so that h[j] and h[j+3] are antipodal
import math
ang = sorted(ring, key=lambda i: math.atan2(float(Sa[i].y), float(Sa[i].x)))
H = ang
print(f"Sa: {n} pts, {len(E)} edges; hexagon {H}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for j in range(3):
    a, c = H[j], H[j + 3]
    assert Sa[a].x == -Sa[c].x and Sa[a].y == -Sa[c].y, "not antipodal"
print("   antipodal pairs:", [(H[j], H[j + 3]) for j in range(3)], flush=True)


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


def survivors(k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    out = []
    for part in partitions(list(range(6))):
        if len(part) > k:
            continue
        # colour symmetry: name the blocks 0, 1, 2, ... in order
        lits = []
        for bi, blk in enumerate(part):
            for e in blk:
                lits.append(1 + H[e] * k + bi)
        if sv.solve(assumptions=lits):
            out.append(tuple(tuple(sorted(x)) for x in part))
    sv.delete()
    return out


for k in (4, 5):
    surv = survivors(k)
    total = sum(1 for p in partitions(list(range(6))) if len(p) <= k)
    anti = [(0, 3), (1, 4), (2, 5)]
    with_mono = [p for p in surv
                 if any(any(x in blk and y in blk for blk in p)
                        for x, y in anti)]
    print(f"\nk={k}: {len(surv)} of {total} hexagon patterns survive"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"   of those, {len(with_mono)} have an antipodal pair together "
          f"-- so the weak property holds iff that is all of them: "
          f"{len(with_mono) == len(surv)}", flush=True)
    byblocks = defaultdict(int)
    for p in surv:
        byblocks[len(p)] += 1
    print(f"   by number of colours used on the hexagon: {dict(sorted(byblocks.items()))}",
          flush=True)
    for p in surv[:12]:
        print("      ", p, flush=True)
    if len(surv) > 12:
        print(f"       ... and {len(surv)-12} more", flush=True)
print("DONE", flush=True)
