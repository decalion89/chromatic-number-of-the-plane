"""Enrich the AUXILIARIES, which is what actually constrains a circle.

Circle points are constrained two ways: by edges among themselves, which are
fixed at 60 degrees apart, and by shared neighbours outside.  A point one away
from two circle points a and b is the pivot reflected across the chord, namely
a + b - p, so the auxiliaries of a 60-point circle are among C(60,2) = 1770
candidates.

Adding the missing ones is the enrichment that could matter: unlike the
pendant circle points, every auxiliary has TWO contacts by construction, so
each one is a genuine constraint linking a pair of circle points.

Question: how many does G already have, how many can be added, and does the
pressure at five colours move once they are in?
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure

g = degrey.build_G()
pts = list(g.vertices)
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
p = pts[p_idx]
circle = sorted(g.adj[p_idx])
print(f"G: {g.n} vertices; pivot {p_idx}, circle of {len(circle)}", flush=True)

have = set(pts)
new = {}
for i, j in itertools.combinations(circle, 2):
    a, b = pts[i], pts[j]
    q = Point(a.x + b.x - p.x, a.y + b.y - p.y)
    if q == p or q in have:
        continue
    new.setdefault(q, set()).update((i, j))
print(f"  {len(new)} auxiliary points missing from G "
      f"(of {len(circle)*(len(circle)-1)//2} chord candidates)", flush=True)
if new:
    counts = sorted((len(v) for v in new.values()), reverse=True)
    print(f"  circle contacts each: {counts[0]} down to {counts[-1]}",
          flush=True)

order = sorted(new, key=lambda q: -len(new[q]))
for keep in (0, 50, 200, len(order)):
    w = build_graph(pts + order[:keep])
    pi = w.index_of(p)
    t = time.time()
    pr = pressure(ColourRelations(w, 5), pi)
    print(f"  +{keep:5} auxiliaries: n={w.n:5} m={w.m:6}, "
          f"pressure at k=5 = {pr}"
          + ("   *** ABOVE 2 ***" if pr > 2 else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
    if keep == len(order):
        break
