"""mu on the 1/sqrt3 ring: the first place it can read above two.

On a circle of radius r, two points are a unit apart at the angle theta with
cos theta = 1 - 1/(2r^2).  A cycle inside the ring needs theta to be a rational
multiple of pi, and cos theta is rational, so Niven's theorem leaves four
radii:

    r = 1         theta =  60 deg   even cycles only   <- the unit circle
    r = 1/sqrt2   theta =  90 deg   even cycles only
    r = 1/sqrt3   theta = 120 deg   TRIANGLES
    r = 1/2       theta = 180 deg   degenerate

Exactly one radius in the plane puts an odd cycle on a ring, and it is
1/sqrt3, the circumradius of a unit triangle.  That is why the number has
turned up in every disjunction this project has built, and it is the reason
the unit circle is hopeless: chi = 2 there, always, in every graph, so a
2-colouring of N(p) exists locally no matter what.

On the 1/sqrt3 ring the local floor is THREE.  So measure mu there:

    mu_5(ring) = min over proper 5-colourings of  |c(ring of radius 1/sqrt3)|

Every mu this project has ever measured read 2.  If this reads 3 it is the
first number above the floor, and the distance to five is two rungs instead of
three.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver

ROOT = HN_DIR
t0 = time.time()
NAMES = sys.argv[1:] or ["five_247_c.json"]
for NAME in NAMES:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    third = F.rational(Fr(1, 3)); one = F.rational(Fr(1))
    ring = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
            if abs(dd - 1/3) < 1e-9 and (g.vertices[i]-g.vertices[j]).norm2() == third:
                ring[i].append(j); ring[j].append(i)
    print(f"\n{NAME} n={n}   [{time.time()-t0:.0f}s]", flush=True)
    ms = MuSolver(g, 5, budget=None)
    print(f"  base colouring {ms.colourable}   [{time.time()-t0:.0f}s]", flush=True)
    got = Counter(); best = None
    order = sorted(ring, key=lambda h: -len(ring[h]))[:120]
    for h in order:
        vs = ring[h]
        tri = sum(1 for a in range(len(vs)) for b in range(a + 1, len(vs))
                  if (g.vertices[vs[a]] - g.vertices[vs[b]]).norm2() == one)
        v = ms.mu(vs)
        got[(v, tri > 0)] += 1
        if best is None or (v or 0) > best[0]:
            best = (v, h, len(vs), tri)
            print(f"    mu = {v} at h={h}, ring {len(vs)} points, "
                  f"{tri} unit edges inside   [{time.time()-t0:.0f}s]", flush=True)
    print(f"  mu on the 1/sqrt3 ring, (value, ring has a triangle): "
          f"{dict(sorted(got.items(), key=str))}   [{time.time()-t0:.0f}s]", flush=True)
    ms.close()
