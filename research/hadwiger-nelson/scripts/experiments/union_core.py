"""The union, measured properly at last.

The requirement for a core at five colours is a 5-chromatic host with the
pivot non-critical, overlapping a pressure gadget. W = G union (G + t) with t
one of the 1/3-vectors meets all three: it contains G so it is 5-chromatic;
removing any vertex leaves a whole copy, so nothing is critical; the
translation folds 224 vertices onto vertices, so the copies are genuinely
coupled; and the three-hexagon pressure gadget is already inside G, since the
hexagon offsets read off the witness are exactly 0, theta/2, theta and their
reflections.

Earlier passes measured this badly. Shrinking from the full target set meant a
first query equivalent to re-proving chi(G) >= 5 with a slow solver; growing
greedily stalled on a plateau. Counterexample construction plus minimisation
costs one solve per point added and at most k per point removed, which is what
this uses.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.forced import (ColourRelations, cegar_core, is_core, minimise_core,
                       pressure)
from hn.geometry import Point
from hn.graph import build_graph

SC = "/tmp/hn"
g = build_G()
mult, (i, j) = pickle.load(open(f"{SC}/topvecs.pkl", "rb"))[0]
dx = g.vertices[j].x - g.vertices[i].x
dy = g.vertices[j].y - g.vertices[i].y
shifted = [Point(v.x + dx, v.y + dy) for v in g.vertices]
pts = list(dict.fromkeys(list(g.vertices) + shifted))
w = build_graph(pts)
gset, sset = set(g.vertices), set(shifted)
overlap = [n for n, v in enumerate(w.vertices) if v in gset and v in sset]
rel = ColourRelations(w, 5)
print(f"W = {w}  5-colourable {rel.colourable}, overlap {len(overlap)}",
      flush=True)

t0, best = time.time(), None
order = sorted(overlap, key=lambda v: -len(w.adj[v])) + \
    sorted(set(range(w.n)) - set(overlap), key=lambda v: -len(w.adj[v]))
for n, bp in enumerate(order):
    c0 = rel.calls
    T, ok = cegar_core(rel, bp, limit=500)
    if not ok:
        if n % 10 == 9:
            print(f"    ... {n+1} pivots, none closed within 500 rounds, "
                  f"best {best[0] if best else None}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
        continue
    raw = len(T)
    T = minimise_core(rel, bp, T)
    if best is None or len(T) < best[0]:
        best = (len(T), bp, [w.vertices[t] for t in T])
        pv = w.vertices[bp]
        d2 = collections.Counter(str(pv.dist2(w.vertices[t])) for t in T)
        print(f"  pivot {bp} deg {len(w.adj[bp])} pressure {pressure(rel, bp)}"
              f": core {raw} -> {len(T)} in {rel.calls-c0} solves, d^2 "
              f"{dict(list(d2.items())[:4])}  [{time.time()-t0:.0f}s]",
              flush=True)
        pickle.dump((bp, best[2]), open(f"{SC}/union_core.pkl", "wb"))
        if len(T) <= 3:
            print(f"  *** CORE OF {len(T)} AT FIVE COLOURS -- BLOCKABLE ***",
                  flush=True)
            break
print(f"best union core: {best[0] if best else None}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
rel.close()
