"""Grow cores from the neighbourhood outwards, at five colours.

T is a core of p exactly when no colour is left unused on N(p) union T -- if
one were, recolouring p to it stays proper and puts c(p) outside c(T). So the
core condition is min |c(N(p) union T)| = k, a pressure measurement, and the
objective is monotone: enlarging T can only raise the minimum. Every previous
core search here shrank a huge target set by repeated separation queries, all
or nothing with no gradient. This grows instead, one solve per candidate,
because a candidate helps exactly when it breaks the current squeeze.

W = G union (G + t) is not vertex-critical, so a core exists. Pressure is 2,
so no core is smaller than three. Three is what the misalignment block closes.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.forced import ColourRelations, cegar_core, is_core, pressure
from hn.geometry import Point
from hn.graph import build_graph

SC = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
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
print(f"W = {w}  5-colourable {rel.colourable}  overlap {len(overlap)}",
      flush=True)

t0, hist, best = time.time(), collections.Counter(), None
order = overlap + sorted(set(range(w.n)) - set(overlap),
                         key=lambda v: -len(w.adj[v]))
for n, bp in enumerate(order):
    T, ok = cegar_core(rel, bp, limit=25)
    if ok:
        hist[len(T)] += 1
        pv = w.vertices[bp]
        d2 = [str(pv.dist2(w.vertices[t])) for t in T]
        if best is None or len(T) < best[0]:
            best = (len(T), bp, T)
            print(f"  *** pivot {bp} deg {len(w.adj[bp])}: CORE OF {len(T)} "
                  f"at d^2 {d2}  [{time.time()-t0:.0f}s]", flush=True)
            pickle.dump((bp, T), open(f"{SC}/grown_{len(T)}_{bp}.pkl", "wb"))
    else:
        hist["no core in 25"] += 1
    if n % 5 == 4:
        print(f"    ... {n+1}/{w.n} pivots, {dict(hist)}, best "
              f"{best[0] if best else None}, {rel.calls} solves  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"cores: {dict(hist)}  best {best[0] if best else None}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
rel.close()
