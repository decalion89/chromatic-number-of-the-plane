"""The ring of radius 1/sqrt3, where the local floor is THREE and not two.

Everything that made the single-point attack hopeless was about the UNIT
circle.  Two points of a radius-1 circle are a unit apart exactly at 60
degrees, so N(p) is a union of paths and hexagons, every cycle is even, and a
2-colouring is just the angle parity.  chi(N(p)) = 2, always, in every graph.

That is a fact about the radius, not about neighbourhoods.  On a circle of
radius r, two points are a unit apart at the angle theta with
cos theta = 1 - 1/(2r^2), and a cycle inside the ring needs theta to be a
rational multiple of pi, so by Niven's theorem cos theta is 0 or +-1/2 or +-1:

    cos theta =  1/2   theta =  60 deg   r = 1        even cycles only
    cos theta =  0     theta =  90 deg   r = 1/sqrt2  even cycles only
    cos theta = -1/2   theta = 120 deg   r = 1/sqrt3  TRIANGLES
    cos theta = -1     theta = 180 deg   r = 1/2      degenerate

One radius in the plane puts an ODD cycle on a ring, and it is 1/sqrt3 -- the
circumradius of a unit triangle, which is why the number has turned up in every
disjunction this project has built.  On that ring the local floor is three
colours, not two, so the gap to five is two rungs instead of three.

What the ring buys is the hub statement, which is the consumable shape:

    forbid c(h) = c(v) for every v at distance 1/sqrt3 from h.
    UNSAT means every 5-colouring puts h's colour somewhere on that ring,

and then rotations about h by 120 degrees carry the ring to itself, every copy
names a partner on the same ring, and two partners 120 degrees apart are a unit
apart -- a contradiction.  SAT is the cheap direction, so the scan is fast and
only a win is slow.
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
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAMES = sys.argv[1:] or ["five_247_c.json", "five_247.json", "five_tuned_1_1.json"]
for NAME in NAMES:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    E = list(g.edges())
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    third = F.rational(Fr(1, 3))
    ring = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
            if abs(dd - 1/3) < 1e-9 and (g.vertices[i]-g.vertices[j]).norm2() == third:
                ring[i].append(j); ring[j].append(i)
    sizes = Counter(len(v) for v in ring.values())
    print(f"\n{NAME} n={n}: {len(ring)} vertices carry a 1/sqrt3 ring; "
          f"sizes {dict(sorted(sizes.items()))}   [{time.time()-t0:.0f}s]", flush=True)
    X = lambda v, c: 1 + v * K + c
    base = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in E:
        for c in range(K):
            base.append([-X(x, c), -X(y, c)])
    # how many colours does the ring itself need, from its own structure?
    tri = 0
    for h, vs in ring.items():
        for a in range(len(vs)):
            for b in range(a + 1, len(vs)):
                u, v = vs[a], vs[b]
                if (g.vertices[u]-g.vertices[v]).norm2() == F.rational(Fr(1)):
                    tri += 1
    print(f"  unit edges inside the rings (triangle sides): {tri}", flush=True)
    best = None; hits = []
    order = sorted(ring, key=lambda h: -len(ring[h]))
    for h in order:
        cnf = list(base)
        for v in ring[h]:
            for c in range(K):
                cnf.append([-X(h, c), -X(v, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(); s.delete()
        if not ok:
            print(f"  *** RING HUB *** h={h} ring={len(ring[h])}   "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            hits.append(h)
            break
    if not hits:
        print(f"  no ring hub: every vertex can avoid its 1/sqrt3 ring   "
              f"[{time.time()-t0:.0f}s]", flush=True)
    else:
        json.dump({"graph": NAME, "hubs": hits},
                  open(f"{ROOT}/data/ring3_{NAME}", "w"))
