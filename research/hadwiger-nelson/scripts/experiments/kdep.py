import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, collections, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_S, build_Sa, build_Y
from hn.forced import ColourRelations, min_colours_on
from hn.graph import build_graph
t0 = time.time()
for name, b in (("S", build_S), ("Sa", build_Sa), ("Y", build_Y)):
    g = b()
    g = g if hasattr(g, "vertices") else build_graph(g)
    for k in (4, 5, 6):
        rel = ColourRelations(g, k)
        if not rel.colourable:
            rel.close(); continue
        h = collections.Counter(min_colours_on(rel, sorted(g.adj[v]))
                                for v in range(g.n) if g.adj[v])
        rel.close()
        print(f"{name:>3} n={g.n:4d} k={k}: max pressure {max(h)} -> smallest "
              f"core {k-max(h)}   {dict(sorted(h.items()))}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
