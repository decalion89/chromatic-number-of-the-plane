"""Where on the curve does the gadget density cross, and what does it demand?

Rigidity is not density, not saturation and not redundancy; what is left is the
NUMBER of independent critical subgraphs per point, each of which imposes one
constraint on a colouring whatever its size.  Sa carries 0.57 Moser spindles
per point and its census is ten; G carries one 5-critical subgraph in 1581
points and its census is the ceiling.

Two things are still unmeasured.  Where along the peel the gadget density
crosses -- the threshold, rather than Sa's own value.  And what that threshold
demands of a construction one level up, which is a matter of arithmetic once
the threshold is known: a union of copies has gadget density exactly
1 / (new points per copy), an invariant of the OVERLAP and not of the number of
copies, so the demand is a bound on how much of each gadget may be new.
"""
import sys, time, math, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
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
rng.shuffle(order)                       # the census curve's own order


def spindles(P):
    bb = IntBasis.covering(P)
    rr = bb.rows(P)
    assert bb.overflow_headroom(rr) < 1.0
    n = len(P)
    dd, DD = bb.dim, bb.D * bb.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(bb, rr)))
    nb = [set() for _ in range(n)]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
    total = 0
    for i in range(n):
        dv = rr - rr[i]
        s = bb._field_square(dv[:, :dd]) + bb._field_square(dv[:, dd:])
        g = s[:, 0] == 3 * DD
        for m in range(1, dd):
            g &= s[:, m] == 0
        far = [int(j) for j in np.nonzero(g)[0]]
        for a in range(len(far) - 1):
            for c in range(a + 1, len(far)):
                x, y = far[a], far[c]
                if y not in nb[x]:
                    continue
                sx, sy = sorted(nb[i] & nb[x]), sorted(nb[i] & nb[y])
                if len(sx) < 2 or len(sy) < 2:
                    continue
                if any(len({i, x, y, sx[p], sx[q], sy[u], sy[v]}) == 7
                       for p in range(len(sx) - 1)
                       for q in range(p + 1, len(sx))
                       for u in range(len(sy) - 1)
                       for v in range(u + 1, len(sy))):
                    total += 1
    return n, total


CENSUS = {107: 715, 157: 715, 207: 715, 257: 715, 307: 715, 347: 577, 397: 10}
print("points  spindles  per point   census at four  [threshold marked]",
      flush=True)
prev = None
for size in sorted(CENSUS):
    P = list(PIN) + order[:size - len(PIN)]
    n, sp = spindles(P)
    mark = ""
    if prev is not None and CENSUS[size] < 715 <= prev:
        mark = "   <- the census leaves the ceiling here"
    prev = CENSUS[size]
    print(f"{n:6d}  {sp:8d}  {sp/n:9.3f}   {CENSUS[size]:6d}{mark}"
          f"  [{time.time()-t0:.0f}s]", flush=True)

print("\nWhat a union of copies can reach:", flush=True)
print("  gadget density of a union = 1 / (new points per copy), which is an",
      flush=True)
print("  invariant of the overlap and NOT of the number of copies.", flush=True)
for label, newpts, size in (("Sa's spindles", 397 / 228, 7),
                            ("G translated (measured)", 1336, 1581),
                            ("a 509-gadget at Sa's overlap fraction",
                             509 * (397 / 228) / 7, 509)):
    print(f"  {label:42s}: {newpts:8.2f} new points per copy -> density "
          f"{1/newpts:.5f}", flush=True)
need = 397 / 228
print(f"\n  To reach Sa's 0.57 with a 509-vertex gadget each copy may cost at "
      f"most {need:.2f} new points,", flush=True)
print(f"  i.e. {100*(1-need/509):.2f} per cent of every gadget must already "
      f"be there.", flush=True)
print("DONE", flush=True)
