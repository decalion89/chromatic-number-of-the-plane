"""Enrich a pivot's circle: the one lever on rho <= deg(p) + 3 not yet pulled.

A core of three at p needs N(p) u T forcing, and supersets of forcing sets
force -- so a BIGGER circle makes the condition easier, not harder.  de Grey's
G has maximum degree 60, but its edge module carries about 134 unit vectors,
so the pivot's circle can in principle hold twice what it does.

The points to add are p + u for each unit step u of the module.  Most are new.
What matters is whether they are CONSTRAINED -- a new circle point adjacent to
nothing else is free, takes any colour, and raises the pressure by nothing.
So each candidate is scored by how many existing vertices it is one away from,
and only the constrained ones are worth adding.

Then the question the whole thing turns on: does the pressure at p rise above
2 once the circle is enriched?
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure

g = degrey.build_G()
pts = list(g.vertices)
one = pts[0].x.field.one()
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
p = pts[p_idx]
print(f"G: {g.n} vertices; pivot {p_idx} has degree {deg[p_idx]}", flush=True)

# the module's unit steps, as differences that occur in the graph
steps = set()
for a, b in g.edges():
    d = (pts[b].x - pts[a].x, pts[b].y - pts[a].y)
    steps.add(d)
    steps.add((-d[0], -d[1]))
print(f"  {len(steps)} unit steps in the edge module", flush=True)

index = {q: i for i, q in enumerate(pts)}
cands = []
for dx, dy in steps:
    q = Point(p.x + dx, p.y + dy)
    if q in index:
        continue
    # how many existing vertices is it one away from?
    touch = sum(1 for r in pts if q.is_unit_apart(r))
    cands.append((touch, q))
cands.sort(key=lambda t: -t[0])
print(f"  {len(cands)} candidate circle points not already present; "
      f"contacts {cands[0][0] if cands else 0} down to "
      f"{cands[-1][0] if cands else 0}", flush=True)

for keep in (0, 10, 30, len(cands)):
    sel = [q for _, q in cands[:keep]]
    w = build_graph(pts + sel)
    pi = w.index_of(p)
    t = time.time()
    pr = pressure(ColourRelations(w, 5), pi)
    print(f"  +{keep:3} circle points: n={w.n:5} m={w.m:6}, "
          f"pivot degree {len(w.adj[pi]):3}, pressure at k=5 = {pr}"
          + ("   *** ABOVE 2 ***" if pr > 2 else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
