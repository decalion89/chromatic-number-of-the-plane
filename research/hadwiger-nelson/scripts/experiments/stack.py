"""Stack copies and watch pressure, which is monotone.

Sa reaches pressure 3 at four colours -- k-1, the classical spindle regime --
and it does so with NO forced pair anywhere on the pivot's circle: the
forced-different graph there is exactly the 30 unit-distance edges, every
degree 2, no odd cycle, no forced non-edge at all. The ambient graph kills
every 2-colouring of the circle jointly, without pinning any single pair.

That matters because it removes the object this package kept failing to find.
Pressure 3 at five colours needs no rainbow and no forced non-edge; it needs
enough ambient rigidity, and pressure measures exactly that, in at most k
solves, and never decreases when points are added. So stack copies and watch.

One translated copy of G gave a flat 2 across nineteen different translations.
This adds several at once, all at |d|^2 = 1/3, the most frequent difference
vector of G and the one carrying 224 vertices onto vertices.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.forced import ColourRelations, min_colours_on
from hn.geometry import Point
from hn.graph import build_graph

SC = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
g = build_G()
top = pickle.load(open(f"{SC}/topvecs.pkl", "rb"))
vecs, seen = [], set()
for mult, (i, j) in top:
    d = (g.vertices[j].x - g.vertices[i].x, g.vertices[j].y - g.vertices[i].y)
    if d not in seen:
        seen.add(d)
        vecs.append(d)
print(f"G = {g}   {len(vecs)} distinct translations available", flush=True)

t0 = time.time()
cur = list(g.vertices)
for copies in range(1, 6):
    if copies > 1:
        dx, dy = vecs[copies - 2]
        cur = list(dict.fromkeys(cur + [Point(q.x + dx, q.y + dy)
                                        for q in cur[:len(g.vertices)]]))
    w = build_graph(cur)
    rel = ColourRelations(w, 5)
    if not rel.colourable:
        print(f"  *** {copies} copies: {w} is NOT 5-COLOURABLE -- "
              f"chi(R^2) >= 6 ***", flush=True)
        rel.close(); break
    hist = collections.Counter()
    order = sorted(range(w.n), key=lambda v: -len(w.adj[v]))
    for v in order:
        hist[min_colours_on(rel, sorted(w.adj[v]))] += 1
    rel.close()
    mx = max(hist)
    print(f"  {copies} copies: {w}  pressures "
          f"{dict(sorted(hist.items()))}  max {mx}  -> smallest core "
          f"{5 - mx}  [{time.time()-t0:.0f}s]", flush=True)
    if mx >= 3:
        print(f"  *** PRESSURE {mx} -- a core of {5-mx} is now possible ***",
              flush=True)
