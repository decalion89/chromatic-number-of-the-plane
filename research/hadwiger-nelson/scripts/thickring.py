"""Thicken de Grey's own ring and ask his own question at five colours.

Rotations about the origin COMMUTE with the sixty-degree rotation, so every
image rho^t(Sa) is dihedrally symmetric exactly as Sa is, and so is any union
of them.  That makes one family the natural continuation of the chain:

    U_m  =  rho^t(Sa),  t = 0..m,  rho = the D = 4 bite

U_0 is Sa, whose D = 4 ring is a hexagon -- six points, three antipodal pairs
-- and whose weak property at FOUR colours is that one of those three pairs is
monochromatic.  U_1 is de Grey's bite, Y, give or take two points.  Beyond
that the ring simply gets thicker: 6(m+1) points and 3(m+1) antipodal pairs,
all at squared distance 16, all sharing one spindle rotation, all about the
one centre.

So the question that broke everywhere else becomes a question about m.  How
thick must the ring be before no five-colouring can keep every antipodal pair
bichromatic?  That is the gateway, one SAT call per level, and a positive
answer is de Grey's chain running at five: weak property on a doubly usable
ring, then the bite, then the spindle.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import doubly_usable_ring
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 14
t0 = time.time()
Sa = build_Sa(K)
ZERO = Point(K.zero(), K.zero())
assert ZERO in Sa, "the origin must be in Sa for its ring to be the hub's"
rho = rotation_joining(4, K)
print(f"Sa: {len(Sa)} points; rho = D=4 bite, doubly usable="
      f"{doubly_usable_ring(Fr(4))}  [{time.time()-t0:.0f}s]", flush=True)

U, seen, cur = list(Sa), set(Sa), list(Sa)
for m in range(0, MAX + 1):
    if m:
        cur = [rho(p) for p in cur]
        fresh = [p for p in cur if p not in seen]
        seen.update(fresh)
        U.extend(fresh)
    b = IntBasis.covering(U)
    r = b.rows(U)
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(U)
    zi = U.index(ZERO)
    d = r - r[zi]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    ring4 = [int(o) for o in np.nonzero(ok)[0]
             if Fr(int(sq[o, 0]), d2) == 4]
    key = {tuple(r[i]): i for i in range(n)}
    rs = set(ring4)
    pairs = []
    for j in ring4:
        anti = key.get(tuple(2 * r[zi] - r[j]))
        if anti is not None and anti in rs and anti > j:
            pairs.append((j, anti))
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    base = sv.solve()
    gate = None
    if base:
        extra = []
        for a, c in pairs:
            for col in range(k):
                extra.append([-(1 + a * k + col), -(1 + c * k + col)])
        sv.append_formula(extra)
        gate = not sv.solve()
    sv.delete()
    tag = "COLOURS" if base else f"*** NOT {k}-COLOURABLE ***"
    gtag = ("*** GATEWAY HOLDS ***" if gate else
            "gateway fails" if gate is not None else "-")
    print(f"m={m:2d}: {n:5d} pts, {len(E):6d} edges, {tag}, D=4 ring "
          f"{len(ring4):3d} pts / {len(pairs):3d} antipodal pairs, {gtag}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not base or gate:
        print("  *** this is the level that closes ***", flush=True)
        break
print("DONE", flush=True)
