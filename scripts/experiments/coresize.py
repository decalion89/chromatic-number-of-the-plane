"""How big is the core, when there is one?

Counterexample-guided construction adds one vertex per solve and terminates
whenever the graph is not k-vertex-critical, so running it without a tight
limit measures the core rather than guessing at it. Twenty-five rounds was
not enough on any pivot of the union; this finds out what is.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.forced import ColourRelations, cegar_core, is_core
from hn.geometry import Point
from hn.graph import build_graph

SC = "/tmp/hn"
g = build_G()
mult, (i, j) = pickle.load(open(f"{SC}/topvecs.pkl", "rb"))[0]
dx, dy = (g.vertices[j].x - g.vertices[i].x, g.vertices[j].y - g.vertices[i].y)
pts = list(dict.fromkeys(list(g.vertices)
                         + [Point(v.x + dx, v.y + dy) for v in g.vertices]))
w = build_graph(pts)
rel = ColourRelations(w, 5)
print(f"W = {w}", flush=True)
t0 = time.time()
for bp in sorted(range(w.n), key=lambda v: -len(w.adj[v]))[:6]:
    c0 = rel.calls
    T, ok = cegar_core(rel, bp, limit=900)
    pv = w.vertices[bp]
    if ok:
        d2 = collections.Counter(str(pv.dist2(w.vertices[t])) for t in T)
        print(f"  pivot {bp} deg {len(w.adj[bp])}: CORE OF {len(T)} in "
              f"{rel.calls-c0} solves, verified {is_core(rel, bp, T)}; "
              f"circles {dict(list(d2.items())[:5])}  [{time.time()-t0:.0f}s]",
              flush=True)
        pickle.dump((bp, T), open(f"{SC}/core_{bp}.pkl", "wb"))
    else:
        print(f"  pivot {bp} deg {len(w.adj[bp])}: no core within 900 "
              f"({rel.calls-c0} solves)  [{time.time()-t0:.0f}s]", flush=True)
rel.close()
