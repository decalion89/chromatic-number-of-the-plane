"""One hexagon at the centre, everything else as dense as it will go.

The ball of radius two over all 62 unit directions the carrier realises gives h
a neighbourhood of 62 points -- and by the grading that is about eleven
components, so 10 * 2^11 = twenty thousand two-colourings to refute.  Richness
in N(h) is the wrong axis, and this is the third time that has shown up.

The burden is minimised at the other extreme: N(h) a SINGLE hexagon, one
component, ten colour pairs times two parities -- twenty patterns, and by colour
symmetry really one question.  Then rung one says exactly

    no 5-colouring makes both inscribed sqrt3-triangles monochromatic.

And a graph is a point set we choose, so we can have both: build the ball over
every direction, then delete the points at distance 1 from h except one hexagon.
The hexagon's own vertices keep their full neighbourhoods, the second shell
keeps everything, and only the thing that was working against us is gone.

That is the sharpest configuration available -- minimum burden at h, maximum
ambient constraint on the triangles that carry it.
"""
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
R = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 9000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(P)
vecs = {}
for x, y in g0.edges():
    for a, b in ((x, y), (y, x)):
        v = Point(g0.vertices[b].x - g0.vertices[a].x,
                  g0.vertices[b].y - g0.vertices[a].y)
        vecs[(round(v.fx, 9), round(v.fy, 9))] = v
V = list(vecs.values())
O = Point(F.rational(0), F.rational(0))
pts = {(0.0, 0.0): O}
q = deque([O])
while q and len(pts) < CAP:
    p = q.popleft()
    for v in V:
        r = Point(p.x + v.x, p.y + v.y)
        if r.fx * r.fx + r.fy * r.fy > R * R + 1e-9: continue
        k = (round(r.fx, 9), round(r.fy, 9))
        if k not in pts:
            pts[k] = r; q.append(r)
one = F.rational(1)
# keep, of the unit circle about h, exactly one hexagon: a direction and its
# five 60-degree rotations, which stay in the field since sqrt3 is in it
r3 = F.sqrt(3); half = F.rational(Fr(1, 2))
c60, s60 = half, r3 * half
best = None
for v in V:
    ring = []
    x, y = v.x, v.y
    ok = True
    for _ in range(6):
        k = (round(float(x), 9), round(float(y), 9))
        if k not in pts: ok = False; break
        ring.append(k)
        x, y = x * c60 - y * s60, x * s60 + y * c60
    if ok and len(set(ring)) == 6:
        best = ring; break
if best is None:
    print("  no full hexagon in the ball", flush=True); sys.exit(0)
keep = {}
for k, p in pts.items():
    dd = p.fx * p.fx + p.fy * p.fy
    if abs(dd - 1.0) < 1e-9 and k not in best: continue
    keep[k] = p
G = build_graph(list(keep.values())); n = G.n
E = list(G.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
h = next(i for i in range(n) if (G.vertices[i] - O).norm2() == F.rational(0))
nb = sorted(adj[h])
inside = sum(1 for a in nb for b in nb
             if a < b and (G.vertices[a]-G.vertices[b]).norm2() == one)
print(f"{NAME}: {len(V)} directions; ball R={R} -> {len(pts)} points, "
      f"{len(pts)-n} unit-circle points dropped", flush=True)
print(f"  graph: n={n} edges={len(E)} degree {2*len(E)/n:.2f}; "
      f"|N(h)|={len(nb)} with {inside} edges inside it   "
      f"[{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])
s = Solver(name="cd19", bootstrap_with=base)
ok = s.solve(); s.delete()
print(f"  5-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if not ok:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"source": NAME, "radius": R, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/ball2_hit.json", "w"))
    sys.exit(0)
ms = MuSolver(G, 5, budget=None)
r = ms.at_most_two(nb)
print(f"  at_most_two(hexagon) = {r}   "
      f"({'placeable' if r else '*** BLOCKED AT TWO: RUNG ONE ***'})   "
      f"[{time.time()-t0:.0f}s]", flush=True)
if r is False:
    print(f"  *** mu_5 = {ms.mu(nb)} ***", flush=True)
    json.dump({"source": NAME, "radius": R, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/ball2_rung1.json", "w"))
ms.close()
