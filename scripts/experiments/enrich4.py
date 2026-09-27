"""The enriched pivot at FOUR colours, as the control for the five-colour result.

Adding every chord point to de Grey's pivot leaves the pressure at five
colours at 2.  That only means something if the same enrichment DOES buy
something at four -- otherwise the enrichment is simply inert and the
five-colour result says nothing about five colours in particular.

So the same graph, the same pivot, at k = 4.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure
from pysat.solvers import Solver

g = degrey.build_G()
pts = list(g.vertices)
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
p = pts[p_idx]
circle = sorted(g.adj[p_idx])
have, extra = set(pts), []
for i, j in itertools.combinations(circle, 2):
    q = Point(pts[i].x + pts[j].x - p.x, pts[i].y + pts[j].y - p.y)
    if q != p and q not in have:
        have.add(q)
        extra.append(q)

# The G rows are dropped: they ask whether G is 4-colourable, which is the
# hard instance and whose answer is already certified. Sa is the control
# that matters -- it IS 4-colourable, so pressure at four colours is a
# real number there rather than a vacuous one.

# and the control that does work: Sa, which IS 4-colourable, enriched the same way
sa = build_graph(degrey.build_Sa())
spts = list(sa.vertices)
sdeg = sa.degrees()
sp = max(range(sa.n), key=lambda v: sdeg[v])
sc = sorted(sa.adj[sp])
shave, sextra = set(spts), []
for i, j in itertools.combinations(sc, 2):
    q = Point(spts[i].x + spts[j].x - spts[sp].x,
              spts[i].y + spts[j].y - spts[sp].y)
    if q != spts[sp] and q not in shave:
        shave.add(q)
        sextra.append(q)
for label, P in (("Sa", spts), ("Sa + all chord points", spts + sextra)):
    w = build_graph(P)
    pi = w.index_of(spts[sp])
    t = time.time()
    pr4 = pressure(ColourRelations(w, 4), pi)
    pr5 = pressure(ColourRelations(w, 5), pi)
    print(f"{label:22}: n={w.n:5} m={w.m:6}  pressure k=4 is {pr4}, "
          f"k=5 is {pr5}  [{time.time()-t:.0f}s]", flush=True)
