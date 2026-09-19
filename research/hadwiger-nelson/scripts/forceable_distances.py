"""Forceable distances at three colours, on the lattice plus its centroids.

On the triangular lattice alone the answer is clean and classical: the
3-colouring is the residue of the Eisenstein integer modulo (1 - omega), of
norm 3, so two points share a colour exactly when the norm of their
difference is divisible by 3. Measured: d^2 = 1, 4 and 7 are forced to
differ, every pair of them, and nothing else is forced at all.

That set does not contain what the three-leg block wants. Its simplest
blocking configuration puts three legs at radius 1/sqrt(3) and angles 0, 180
and 60 degrees, whose pairwise squared distances are 1, 1/3 and 4/3 -- and
1/3 and 4/3 are not lattice distances. Centroids of the lattice's triangles
sit at exactly 1/sqrt(3) from their three corners, so adding them is what
puts those distances in the graph at all.
"""
import sys, time, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.field import Field
from hn.forced import ColourRelations
from hn.geometry import Point
from hn.graph import build_graph

f = Field((3,))
r3 = f.sqrt(3)
half, sixth, third = (f.rational(Fr(1, 2)), f.rational(Fr(1, 6)),
                      f.rational(Fr(1, 3)))
steps = [Point(f.rational(1), f.zero()),
         Point(half, r3 * half),
         Point(-half, r3 * half)]
lat = {Point(f.zero(), f.zero())}
for _ in range(3):
    new = set()
    for q in lat:
        for s in steps:
            new.add(Point(q.x + s.x, q.y + s.y))
            new.add(Point(q.x - s.x, q.y - s.y))
    lat |= new
# centroids: (u + v + w)/3 of every unit triangle
cent = set()
L = sorted(lat, key=lambda q: float(q.x * q.x + q.y * q.y))
gl = build_graph(L)
one = f.rational(1)
for u in range(gl.n):
    for v in gl.adj[u]:
        for w in gl.adj[u] & gl.adj[v]:
            if u < v < w:
                cent.add(Point((gl.vertices[u].x + gl.vertices[v].x
                                + gl.vertices[w].x) * third,
                               (gl.vertices[u].y + gl.vertices[v].y
                                + gl.vertices[w].y) * third))
pts = sorted(lat | cent, key=lambda q: float(q.x * q.x + q.y * q.y))
g = build_graph(pts)
print(f"lattice {len(lat)} + centroids {len(cent)} -> {g}", flush=True)

rel = ColourRelations(g, 3)
print(f"3-colourable: {rel.colourable}", flush=True)
if not rel.colourable:
    sys.exit(0)
t0, forced, free, checked = time.time(), collections.Counter(), collections.Counter(), 0
for a in range(g.n):
    for b in range(a + 1, g.n):
        d2 = g.vertices[a].dist2(g.vertices[b])
        if float(d2) > 5.0:
            continue
        checked += 1
        key = str(d2)
        (forced if rel.different(a, b) else free)[key] += 1
    if a % 50 == 49:
        print(f"    ... {a+1}/{g.n}, {checked} pairs, {len(forced)} forceable"
              f"  [{time.time()-t0:.0f}s]", flush=True)
rel.close()
print(f"\n{checked} pairs.  forced-to-differ squared distances:", flush=True)
for key in sorted(forced, key=lambda s: -forced[s]):
    print(f"   d^2 = {key:>26}  forced {forced[key]:5d}  free {free.get(key,0):5d}",
          flush=True)
print(f"\nnever forced (sample): "
      f"{sorted(set(free) - set(forced), key=lambda s: free[s], reverse=True)[:10]}",
      flush=True)
