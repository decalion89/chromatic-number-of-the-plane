"""Symmetrise G, the way de Grey symmetrised S.

The cap test needs two things at once and no object in this work has both.

It needs a TIGHT carrier -- chromatic number equal to the number of colours
asked about -- or a spare colour makes every ring showable.  And it needs
RINGS: a centre with a full orbit of points at one distance, for the palette
bound to bite on.

Sa has both at four colours: it is 4-chromatic, and it is the closure of the
39 seed points under the order-twelve group of the hexagonal rotation and the
mirror, so every point sits on a complete ring about the origin.  That is not
incidental to de Grey's construction, it is its first step.

At five colours the only tight carrier available is G, and G has no symmetry
at all: it is built by turning Y about the point (-2, 0), so the order-twelve
group about the origin does not act on it, and its 1581 points fall into 1581
orbits of size one.  No symmetry, no rings, no cap.

So take de Grey's first step one level up and close G under the same group.
The result contains G, so it needs at least five colours; and every point of
it sits on a complete ring about the origin by construction.  If it needs six,
that is the answer.  If it needs five, it is the first carrier in this work
that is tight and symmetric at once, and the palette depth measured on its
rings is the real analogue of de Grey's lemma.
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
rot = _rot60(K)
base = build_G(K, as_graph=False)
seen, P = set(), []
for p in base:
    for q0 in (p, Point(p.x, -p.y)):
        q = q0
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                P.append(q)
            q = rot(q)
print(f"symmetric closure of G: {len(P)} points from {len(base)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
b = IntBasis.covering(P)
r = b.rows(P)
hr = b.overflow_headroom(r)
print(f"overflow headroom {hr:.3f} (must be < 1)  [{time.time()-t0:.0f}s]",
      flush=True)
assert hr < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{n} points, {len(E)} edges ({2*len(E)/n:.2f}/v)"
      f"  [{time.time()-t0:.0f}s]", flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/symG.pkl",
          "wb") as f:
    pickle.dump(P, f)

# Rings about the origin exist by construction; report the biggest.
dm, D2 = b.dim, b.D * b.D
oi = next((i for i, p in enumerate(P)
           if not any(p.x.c) and not any(p.y.c)), 0)
d = r - r[oi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
g = defaultdict(int)
for off in np.nonzero(ok)[0]:
    if off != oi:
        D = Fr(int(sq[off, 0]), D2)
        if 0 < D <= 4:
            g[D] += 1
big = sorted(((c, str(D)) for D, c in g.items()), reverse=True)[:8]
print(f"rings about the origin: {big}  [{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
for col in range(1, k):
    cls.append([-(1 + col)])
print(f"solving at k={k}, {n*k} variables, {len(cls)} clauses"
      f"  [{time.time()-t0:.0f}s]", flush=True)
s = Solver(name="cd15", bootstrap_with=cls)
ok2 = s.solve()
s.delete()
if ok2:
    print(f"{k}-colourable: True  [{time.time()-t0:.0f}s]", flush=True)
else:
    print(f"*** NOT {k}-COLOURABLE -- chi > {k} *** "
          f" [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
