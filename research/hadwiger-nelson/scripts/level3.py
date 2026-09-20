"""One level further out: complete every unit edge into equilateral triangles.

Local enrichment around the pivot is exhausted -- every chord point added,
1680 auxiliaries, pressure still 2 and the confined set still 2-degenerate.
What that cannot see is structure beyond the first two levels.

The cheapest way further out: for every unit-distance pair u, v in the graph,
the two points completing an equilateral triangle on uv.  They exist in de
Grey's field because sqrt 3 does, they are one away from BOTH u and v by
construction, and they are exactly the points a spindle would reach for.

Then the only question that matters: does the pressure at the pivot move off
2 at five colours?
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure

g = degrey.build_G()
pts = list(g.vertices)
fld = pts[0].x.field
HALF = fld.rational(Fraction(1, 2))
RT3 = fld.sqrt(3) * HALF
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
p = pts[p_idx]
circle = sorted(g.adj[p_idx])

# first the chord enrichment, as before
have, extra = set(pts), []
for i, j in itertools.combinations(circle, 2):
    q = Point(pts[i].x + pts[j].x - p.x, pts[i].y + pts[j].y - p.y)
    if q != p and q not in have:
        have.add(q)
        extra.append(q)
w = build_graph(pts + extra)
print(f"chord-enriched: {w.n} vertices, {w.m} edges", flush=True)

# then equilateral completions of every unit edge, near the pivot only --
# the whole graph would be 30000 points and the pivot's pressure cannot see
# structure far away
near = [v for v in range(w.n) if float(w.vertices[v].dist2(p)) <= 9.0]
nearset = set(near)
print(f"  {len(near)} vertices within 3 of the pivot", flush=True)
t = time.time()
new = set()
for a in near:
    ua = w.vertices[a]
    for b in w.adj[a]:
        if b <= a or b not in nearset:
            continue
        ub = w.vertices[b]
        mx, my = (ua.x + ub.x) * HALF, (ua.y + ub.y) * HALF
        dx, dy = ub.x - ua.x, ub.y - ua.y
        for s in (1, -1):
            q = Point(mx - dy * RT3 * fld.rational(s),
                      my + dx * RT3 * fld.rational(s))
            if q not in have:
                new.add(q)
print(f"  {len(new)} equilateral completions found  "
      f"[{time.time()-t:.0f}s]", flush=True)

order = sorted(new, key=lambda q: float(q.dist2(p)))
for keep in (0, 200, 1000, len(order)):
    z = build_graph(pts + extra + order[:keep])
    pi = z.index_of(p)
    t = time.time()
    pr = pressure(ColourRelations(z, 5), pi)
    print(f"  +{keep:5} completions: n={z.n:6} m={z.m:7}, "
          f"pivot degree {len(z.adj[pi]):3}, pressure k=5 = {pr}"
          + ("   *** ABOVE 2 ***" if pr > 2 else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
    if keep == len(order):
        break
