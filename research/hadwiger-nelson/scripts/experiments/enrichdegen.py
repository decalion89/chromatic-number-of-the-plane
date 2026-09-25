"""Does the fully enriched auxiliary graph escape 2-degeneracy?  No.

Adding every chord point a + b - p to de Grey's pivot -- 1380 new auxiliaries,
each with exactly two circle contacts, taking the graph from 1581 to 2961
vertices -- leaves the pressure at five colours unchanged at 2.  That is the
degeneracy result arriving from a completely different direction, so it is
worth confirming in its own terms: is the enriched confined set still
2-degenerate?

If it is, the earlier hexagon-family scans were not measuring a quirk of that
family.  Maximal local enrichment, in the field the construction lives in,
cannot produce the degeneracy 3 that lists of size 3 need.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph


def degeneracy(adj, verts):
    rem = {v: len(adj[v] & verts) for v in verts}
    d = 0
    while rem:
        v = min(rem, key=rem.get)
        d = max(d, rem[v])
        for u in adj[v]:
            if u in rem:
                rem[u] -= 1
        del rem[v]
    return d


g = degrey.build_G()
pts = list(g.vertices)
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
p = pts[p_idx]
circle = sorted(g.adj[p_idx])
have = set(pts)
extra = []
for i, j in itertools.combinations(circle, 2):
    q = Point(pts[i].x + pts[j].x - p.x, pts[i].y + pts[j].y - p.y)
    if q != p and q not in have:
        have.add(q)
        extra.append(q)
w = build_graph(pts + extra)
pi = w.index_of(p)
circ = w.adj[pi]
aux = {v for v in range(w.n)
       if v != pi and v not in circ and len(w.adj[v] & circ) >= 2}
print(f"enriched: {w.n} vertices, {w.m} edges; circle {len(circ)}, "
      f"{len(aux)} auxiliaries", flush=True)
print(f"  auxiliary graph degeneracy: {degeneracy(w.adj, aux)}", flush=True)

# the circle splits into components; enumerate its proper 2-colourings
cadj = {v: w.adj[v] & circ for v in circ}
comps, seen = [], set()
for v in circ:
    if v in seen:
        continue
    stack, part = [v], []
    seen.add(v)
    while stack:
        u = stack.pop()
        part.append(u)
        for x in cadj[u]:
            if x not in seen:
                seen.add(x)
                stack.append(x)
    comps.append(part)
base = []
for part in comps:
    col, stack, ok = {}, [(part[0], 0)], True
    while stack:
        u, c = stack.pop()
        if u in col:
            ok &= col[u] == c
            continue
        col[u] = c
        for x in cadj[u]:
            stack.append((x, 1 - c))
    if not ok:
        print("  a circle component is not bipartite -- pressure > 2 free")
        sys.exit()
    base.append([col, {u: 1 - c for u, c in col.items()}])
print(f"  circle in {len(comps)} components, "
      f"{2 ** len(comps)} orientations", flush=True)

worst = 0
for choice in itertools.product(*base):
    col = {}
    for d in choice:
        col.update(d)
    conf = {v for v in aux if {col[x] for x in w.adj[v] & circ} == {0, 1}}
    worst = max(worst, degeneracy(w.adj, conf))
print(f"  worst confined-set degeneracy over all orientations: {worst}"
      + ("   *** 3 OR MORE ***" if worst >= 3 else "   (still 2 or less)"),
      flush=True)
