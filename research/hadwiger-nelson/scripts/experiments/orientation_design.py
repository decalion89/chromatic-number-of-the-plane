"""Build the gadget the lift needs, instead of hunting for it.

The mechanism is now named: the pivot's circle is five hexagons, five free
orientation bits; each auxiliary point is adjacent to two circle points in
DIFFERENT hexagons, so it is barred from both colours exactly when those two
hexagons are oriented against each other; and the barred points, confined to
k-2 colours, must fail to be colourable. At four colours an odd cycle does it.
At five the confined set must be 4-chromatic.

One thing the mechanism rules out immediately. A point whose two circle
neighbours are ADJACENT -- 60 degrees apart, same hexagon -- is barred
whatever the orientation, and those points are exactly the far apexes of the
unit rhombi, all at distance sqrt(3) from the pivot. On that circle the
chord-1 angle is 2 arcsin(1/(2 sqrt 3)) = 33.56 degrees, which does not close,
so their graph is a union of paths: bipartite, never even an odd cycle. The
orientation-dependence is not an accident of de Grey's construction, it is
forced.

So: take the witness circle, generate every point one away from two circle
points of different hexagons, and for each of the 32 orientations ask what the
barred set looks like. The number that matters is its chromatic number.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.graph import build_graph
from hn.mixed import circle_intersections

R = HN_DIR
pts, doc = load_certificate(f"{R}/certificates/pressure3_witness_47.json")
g = build_graph(pts)
field = pts[0].x.field
one = field.rational(1)
piv = 0
circle = sorted(g.adj[piv])
seen, comps = set(), []
for s in circle:
    if s in seen:
        continue
    comp, st = [], [s]
    seen.add(s)
    while st:
        x = st.pop()
        comp.append(x)
        for y in g.adj[x]:
            if y in set(circle) and y not in seen:
                seen.add(y)
                st.append(y)
    comps.append(sorted(comp))
print(f"circle {len(circle)} in {len(comps)} hexagons", flush=True)

# parity within each hexagon: alternate around the 6-cycle
par = {}
for ci, comp in enumerate(comps):
    start = comp[0]
    par[start] = 0
    front = [start]
    while front:
        x = front.pop()
        for y in g.adj[x]:
            if y in set(comp) and y not in par:
                par[y] = 1 - par[x]
                front.append(y)
comp_of = {v: ci for ci, comp in enumerate(comps) for v in comp}

# every exact point one away from two circle points of different hexagons
t0 = time.time()
cand, seenq = [], set(pts)
for a, b in itertools.combinations(circle, 2):
    if comp_of[a] == comp_of[b]:
        continue
    if float(g.vertices[a].dist2(g.vertices[b])) > 3.999:
        continue
    for q in circle_intersections(g.vertices[a], one, g.vertices[b], one):
        if q in seenq:
            continue
        seenq.add(q)
        cand.append((q, a, b))
print(f"{len(cand)} candidate auxiliary points  [{time.time()-t0:.0f}s]",
      flush=True)

# each candidate is barred exactly when its two circle neighbours differ:
#   colour(v) = par[v] XOR orientation[comp_of[v]]
qs = [q for q, _a, _b in cand]
gq = build_graph(qs)
idx = {q: i for i, q in enumerate(gq.vertices)}
worst = None
for bits in range(1 << len(comps)):
    def col(v):
        return par[v] ^ ((bits >> comp_of[v]) & 1)
    barred = [idx[q] for q, a, b in cand if col(a) != col(b)]
    sub = build_graph([gq.vertices[i] for i in barred])
    chi = 2
    while chi < 6 and not is_k_colorable(sub, chi)[0]:
        chi += 1
    if worst is None or chi < worst[0]:
        worst = (chi, bits, len(barred))
        print(f"  orientation {bits:05b}: {len(barred)} barred, chromatic "
              f"number {chi}  [{time.time()-t0:.0f}s]", flush=True)
print(f"WORST orientation: chromatic number {worst[0]} over {worst[2]} barred "
      f"points.  Need >= 3 for k=4 (have it) and >= 4 for k=5  "
      f"[{time.time()-t0:.0f}s]", flush=True)
