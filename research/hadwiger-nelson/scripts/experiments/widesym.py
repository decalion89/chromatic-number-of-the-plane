"""Close G under the group about a FAR centre, to manufacture a wide ring.

The third level of de Grey's recursion wants a ring of radius 8, and nothing
built here has one: a large ring needs many points at exactly that distance
from a single centre, and a graph of diameter 6 does not spread that far.
Symmetry is what puts many points on one circle -- Sa's six at radius 2, Gp's
twelve at radius 4 are orbits, not accidents -- so the ring has to be made,
not found.

Closing about a far centre makes it.  Twelve images of G with their centres
on a circle of radius rho sit 2*rho*sin(15 degrees) apart, about 0.52*rho,
and each has diameter about 6, so they still overlap as long as rho is under
roughly 11.  At rho = 8 the copies overlap AND every point of G at distance
exactly 8 from the centre generates a twelve-point ring there -- which is
arranged simply by putting the centre 8 away from a point of G.

That is the opposite choice from Gp, where the centre was the pivot two units
out and the copies piled on top of each other.  Far enough to reach, close
enough to touch.

The bite on a radius-8 ring is cos 127/128 with sin sqrt(3*5*17)/128 and
needs sqrt17, but the closure itself does not, so the carrier is built and
measured first.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
OFF = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(8)
t0 = time.time()
base = build_G(K, as_graph=False)
anchor = base[0]
C = Point(anchor.x + K.rational(OFF), anchor.y)
print(f"anchor ({float(anchor.x):.3f},{float(anchor.y):.3f}), centre "
      f"({float(C.x):.3f},{float(C.y):.3f}), offset {OFF}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
rot = _rot60(K).about(C)


def mirror(p):
    return Point(p.x, C.y + (C.y - p.y))


seen, P = set(), []
for p in base:
    for q0 in (p, mirror(p)):
        q = q0
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                P.append(q)
            q = rot(q)
n = len(P)
print(f"closed about the far centre: {n} points from {len(base)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
b = IntBasis.covering(P)
r = b.rows(P)
hr = b.overflow_headroom(r)
print(f"overflow headroom {hr:.3f}  [{time.time()-t0:.0f}s]", flush=True)
assert hr < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
print(f"{n} points, {len(E)} edges ({2*len(E)/n:.2f}/v; G itself has "
      f"{2*7877/1581:.2f}/v, Gp has 10.64/v)  [{time.time()-t0:.0f}s]",
      flush=True)
with open(f"/tmp/hn/wide{OFF}.pkl",
          "wb") as f:
    pickle.dump(P, f)

rc = b.rows([C])[0]
dm, D2 = b.dim, b.D * b.D
d = r - rc
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
g = defaultdict(int)
for off in np.nonzero(ok)[0]:
    g[Fr(int(sq[off, 0]), D2)] += 1
big = sorted(((c, str(Dv)) for Dv, c in g.items()), reverse=True)[:10]
print(f"rings about the far centre: {big}  [{time.time()-t0:.0f}s]",
      flush=True)
tgt = Fr(OFF) ** 2
print(f"the radius-{OFF} ring (D={tgt}) holds {g.get(tgt, 0)} points",
      flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
for col in range(1, k):
    cls.append([-(1 + col)])
print(f"solving at k={k}  [{time.time()-t0:.0f}s]", flush=True)
s = Solver(name="cd15", bootstrap_with=cls)
ok2 = s.solve()
s.delete()
print(f"{k}-colourable: {ok2}" if ok2 else
      f"*** NOT {k}-COLOURABLE -- chi > {k} ***", flush=True)
print(f"  [{time.time()-t0:.0f}s]\nDONE", flush=True)
