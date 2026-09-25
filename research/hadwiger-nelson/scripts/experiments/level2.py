"""A second level: points attached to the auxiliaries, not to the circle.

Everything one away from two circle points is a Minkowski sum u + v, so the
depth-one design is already exhausted -- no more confined points exist without
adding circle points. The list formulation stops there too: a point attached
to AUXILIARIES has no list, because the auxiliaries are not fixed, only
constrained. What it has is ordinary propagation, and only a solver sees it.

So: take the three-hexagon gadget, add every exact point one away from two of
its auxiliaries, and measure the pressure again. Pressure is monotone, so this
can only help, and it is the last lever the geometry offers before enlarging
the circle.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.forced import ColourRelations, min_colours_on
from hn.geometry import Point
from hn.graph import build_graph
from hn.mixed import circle_intersections, three_hexagon_gadget

field, pivot, pts = three_hexagon_gadget()
one = field.rational(1)
g = build_graph(pts)
piv = g.vertices.index(pivot)
circle = set(g.adj[piv])
aux = [v for v in range(g.n) if v != piv and v not in circle]
print(f"base {g}: circle {len(circle)}, auxiliaries {len(aux)}", flush=True)

t0 = time.time()
have = set(pts)
cand = {}
for a, b in itertools.combinations(aux, 2):
    if float(g.vertices[a].dist2(g.vertices[b])) > 3.999:
        continue
    for q in circle_intersections(g.vertices[a], one, g.vertices[b], one):
        if q not in have:
            cand[q] = cand.get(q, 0) + 1
print(f"{len(cand)} level-two candidates  [{time.time()-t0:.0f}s]", flush=True)

ranked = sorted(cand, key=lambda q: -cand[q])
for take in (0, 200, 600, 1200, len(ranked)):
    extra = ranked[:take]
    h = build_graph(pts + extra)
    j = h.vertices.index(pivot)
    rel = ColourRelations(h, 5)
    if not rel.colourable:
        print(f"  +{take}: {h} is NOT 5-COLOURABLE -- chi(R^2) >= 6 ***",
              flush=True)
        rel.close()
        break
    t = min_colours_on(rel, sorted(h.adj[j]))
    rel.close()
    print(f"  +{take} level-two points: {h}  pressure at k=5: {t}"
          f"{'   *** PRESSURE 3 ***' if t >= 3 else ''}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if t >= 3:
        import pickle
        pickle.dump(pts + extra, open("/tmp/claude-0/-home-user-darwin-50/"
                                      "aceaa9ec-f432-5848-a506-39c59179b415/"
                                      "scratchpad/pressure3_k5.pkl", "wb"))
        break
