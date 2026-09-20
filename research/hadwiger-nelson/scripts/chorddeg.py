"""How many of the added chord points are individually inert?

A vertex with two neighbours is colourable whenever k >= 3, so a chord point
whose ONLY neighbours are two circle points constrains nothing at all.  The
pressure control said as much: adding every chord point to Sa leaves the
pressure unchanged at both four and five colours.

So the honest question is not how many were added but how many acquired more
than two neighbours -- those are the only ones that can carry a constraint,
and their mutual adjacencies are what the confined-set degeneracy measures.
"""
import sys, itertools, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph

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
w = build_graph(pts + extra)
wdeg = w.degrees()
idx = {q: i for i, q in enumerate(w.vertices)}
added = [idx[q] for q in extra]
dist = collections.Counter(wdeg[v] for v in added)
print(f"{len(added)} chord points added; degree distribution "
      f"{dict(sorted(dist.items()))}", flush=True)
inert = sum(n for d, n in dist.items() if d <= 2)
print(f"  {inert} of {len(added)} have degree <= 2 and are individually "
      f"inert at any k >= 3  ({100*inert/len(added):.0f} per cent)", flush=True)
pi = idx[p]
circ = w.adj[pi]
useful = [v for v in added if wdeg[v] > 2]
print(f"  {len(useful)} carry more than two neighbours; among themselves they "
      f"span {sum(1 for a in useful for b in w.adj[a] if b in set(useful)) // 2}"
      f" edges", flush=True)
