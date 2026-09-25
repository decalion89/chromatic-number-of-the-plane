"""Close G under the order-twelve group about its own centre, not the origin.

Symmetrising is de Grey's first move and it only works at the right centre.
S straddles the origin, so closing it there makes Sa overlap itself and become
a single dense graph.  G does not: it is Y turned about the point (-2, 0), it
lives entirely in the upper half plane around (-2, 2), and closing it about
the origin produced twelve copies so far apart that the union had exactly
G's own edge density, 9.97 per vertex -- twelve disjoint graphs wearing one
name, 5-colourable for free.

The centre G actually has is the pivot (-2, 0).  It is a vertex of G, 129 of
G's 1581 points sit at rational distance from it in 31 rings, and the twelve
images of G about it are centred two apart on a circle of radius two while
each has radius about 2.8 -- so they overlap heavily instead of standing
apart.  That is the geometry that makes a symmetrisation worth doing.

What comes out contains G, so it needs five colours.  It is closed under the
group, so every point sits on a complete ring about the pivot.  Tight and
symmetric at once, which is de Grey's situation at four colours, one level up.
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
t0 = time.time()
PIV = Point(K.rational(-2), K.zero())
rot = _rot60(K).about(PIV)


def mirror(p):                      # reflection in the line y = 0
    return Point(p.x, -p.y)


base = build_G(K, as_graph=False)
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
print(f"G closed about the pivot: {n} points from {len(base)} "
      f"(12 x {len(base)} = {12*len(base)} before overlap)"
      f"  [{time.time()-t0:.0f}s]", flush=True)
b = IntBasis.covering(P)
r = b.rows(P)
hr = b.overflow_headroom(r)
print(f"overflow headroom {hr:.3f}  [{time.time()-t0:.0f}s]", flush=True)
assert hr < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
print(f"{n} points, {len(E)} edges ({2*len(E)/n:.2f}/v; G itself has "
      f"{2*7877/1581:.2f}/v)  [{time.time()-t0:.0f}s]", flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/pivG.pkl",
          "wb") as f:
    pickle.dump(P, f)

pi = next(i for i, p in enumerate(P) if p == PIV)
dm, D2 = b.dim, b.D * b.D
d = r - r[pi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
g = defaultdict(int)
for off in np.nonzero(ok)[0]:
    if off != pi:
        D = Fr(int(sq[off, 0]), D2)
        if 0 < D <= 4:
            g[D] += 1
print(f"rings about the pivot (D <= 4): "
      f"{sorted(((c, str(D)) for D, c in g.items()), reverse=True)[:10]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
sym = list(cls)
for col in range(1, k):
    sym.append([-(1 + col)])
print(f"solving at k={k}: {n*k} variables, {len(sym)} clauses"
      f"  [{time.time()-t0:.0f}s]", flush=True)
s = Solver(name="cd15", bootstrap_with=sym)
ok2 = s.solve()
s.delete()
print(f"{k}-colourable: {ok2}" if ok2 else
      f"*** NOT {k}-COLOURABLE -- chi > {k} ***", flush=True)
print(f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
