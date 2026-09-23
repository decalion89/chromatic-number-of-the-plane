"""Adjoin the centre of every unit triangle -- the operation that makes 1/sqrt3.

The win condition is proved and narrow: if some vertex h has two points v1, v2 at
distance 1/sqrt3 from it with c(h) = c(v1) or c(h) = c(v2) in every 5-colouring,
then three copies rotated by 120 degrees about h are 6-chromatic.  Two copies
that name the same index land on points of the same 1/sqrt3 circle 120 degrees
apart, which is distance 1, and pigeonhole on two indices across three copies
forces that.  Nothing else about the geometry is used.

The necessary first step is a graph where a vertex cannot avoid its own 1/sqrt3
ring at all, and none of the graphs here can: every vertex escapes, rings of
eighteen points included.  Scanning more of them is pointless -- what is needed
is a graph richer in the relation itself.

1/sqrt3 is the circumradius of a unit triangle, so the operation that
manufactures the relation is to adjoin every unit triangle's CENTRE.  Each new
point arrives at distance exactly 1/sqrt3 from three mutually adjacent old ones,
which is a whole triangle on its ring at once, and the centroid of three field
points stays in the field, so nothing is extended.

Three questions, cheapest first: is the enriched graph still 5-colourable; can
any vertex still avoid its ring; and does forbidding both distances at once now
kill it.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from itertools import combinations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n0 = g.n
adj = defaultdict(set)
for x, y in g.edges():
    adj[x].add(y); adj[y].add(x)
third = F.rational(Fr(1, 3))
tri = []
for u in range(n0):
    for v in sorted(adj[u]):
        if v <= u: continue
        for w in sorted(adj[u] & adj[v]):
            if w > v: tri.append((u, v, w))
print(f"{NAME} n={n0}: {len(tri)} unit triangles   [{time.time()-t0:.0f}s]",
      flush=True)
onethird = F.rational(Fr(1, 3))
pts = {(round(q.fx, 9), round(q.fy, 9)): q for q in g.vertices}
added = 0
for u, v, w in tri:
    c = Point((g.vertices[u].x + g.vertices[v].x + g.vertices[w].x) * onethird,
              (g.vertices[u].y + g.vertices[v].y + g.vertices[w].y) * onethird)
    k = (round(c.fx, 9), round(c.fy, 9))
    if k not in pts:
        pts[k] = c; added += 1
U = build_graph(list(pts.values()))
print(f"  + {added} triangle centres -> n={U.n}, edges={sum(1 for _ in U.edges())}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

n = U.n; E = list(U.edges())
hx = [q.fx for q in U.vertices]; hy = [q.fy for q in U.vertices]
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])
s = Solver(name="cd19", bootstrap_with=base)
ok = s.solve()
print(f"  5-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if not ok:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"source": NAME, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in U.vertices]},
              open(f"{ROOT}/data/centres_hit.json", "w"))
    sys.exit(0)
s.delete()

ring = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if abs(dd - 1/3) < 1e-9 and (U.vertices[i]-U.vertices[j]).norm2() == third:
            ring[i].append(j); ring[j].append(i)
print(f"  {len(ring)} vertices carry a 1/sqrt3 ring; sizes "
      f"{dict(sorted(Counter(len(v) for v in ring.values()).items()))}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
order = sorted(ring, key=lambda h: -len(ring[h]))
for h in order[:60]:
    cnf = list(base)
    for v in ring[h]:
        for c in range(K):
            cnf.append([-X(h, c), -X(v, c)])
    sv = Solver(name="cd19", bootstrap_with=cnf)
    good = sv.solve(); sv.delete()
    if not good:
        print(f"  *** RING HUB h={h}, ring {len(ring[h])} ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        json.dump({"source": NAME, "hub": h, "ring": ring[h]},
                  open(f"{ROOT}/data/ringhub_hit.json", "w"))
        break
else:
    print(f"  no ring hub among the 60 richest   [{time.time()-t0:.0f}s]", flush=True)

pairs = [(i, j) for i in ring for j in ring[i] if j > i]
cnf = list(base)
for i, j in pairs:
    for c in range(K):
        cnf.append([-X(i, c), -X(j, c)])
sv = Solver(name="cd19", bootstrap_with=cnf)
r = sv.solve(); sv.delete()
print(f"  chi(U; 1 and 1/sqrt3) <= 5 with {len(pairs)} extra pairs: {r}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
