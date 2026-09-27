"""Local density instead of global size: the densest ball around the target.

Rung one asks that the hexagon about h never be 2-coloured, which means its two
inscribed sqrt3-triangles are never both monochromatic.  If T is monochromatic
in gamma then N(T) -- everything a unit from T -- loses gamma, and likewise
N(T') loses delta.  So the hypothesis costs whatever fraction of the graph those
two regions cover.

In every graph built here, that fraction is tiny.  The hexagon's six vertices
have degree about thirty, so N(T) u N(T') is a couple of hundred points out of
four thousand -- five per cent of the graph loses one colour each, and the rest
absorbs it without noticing.  Enriching globally makes the denominator grow too,
which is why enrichment never helped.

So invert the construction.  Take the unit vectors the carriers actually
realise, and grow the BALL around h: every point reachable from h by unit steps
within radius R, and nothing else.  Then the graph is exactly the part of the
plane that the hexagon's two triangles can reach, N(T) u N(T') is most of it,
and the 2-colouring hypothesis has nowhere to hide.

Small, dense, and centred on the question -- the opposite of every graph in
data/, which are large, sparse, and centred nowhere.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
R = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 6000
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
print(f"{NAME}: {len(V)} distinct unit directions   [{time.time()-t0:.0f}s]",
      flush=True)
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
G = build_graph(list(pts.values())); n = G.n
E = list(G.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
h = next(i for i in range(n) if (G.vertices[i] - O).norm2() == F.rational(0))
print(f"  ball of radius {R}: n={n} edges={len(E)} degree {2*len(E)/n:.2f}; "
      f"deg(h)={len(adj[h])}   [{time.time()-t0:.0f}s]", flush=True)
one = F.rational(1); three = F.rational(3)
nb = sorted(adj[h])
hexa = [v for v in nb
        if sum(1 for w in nb if (G.vertices[v]-G.vertices[w]).norm2() == one) >= 2]
print(f"  |N(h)|={len(nb)}, of which {len(hexa)} have two neighbours inside it   "
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
              open(f"{ROOT}/data/ball_hit.json", "w"))
    sys.exit(0)
ms = MuSolver(G, 5, budget=None)
r = ms.at_most_two(nb)
print(f"  at_most_two(N(h)) = {r}   "
      f"({'placeable' if r else '*** BLOCKED AT TWO: RUNG ONE ***'})   "
      f"[{time.time()-t0:.0f}s]", flush=True)
if r is False:
    v = ms.mu(nb)
    print(f"  *** mu_5(N(h)) = {v} ***", flush=True)
    json.dump({"source": NAME, "radius": R, "mu": v,
               "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/ball_rung1.json", "w"))
ms.close()
