"""Put the critical structure ON the hexagon, not near it.

The ideal shape is now built and measured: a single hexagon at h, so one
component and twenty patterns; a dense ball around it; and, merged in, a
translate of the vertex-critical 5-chromatic graph so the whole thing refuses
four.  7704 points, degree 9.23, |N(h)| = 6.  And the hexagon still 2-colours.

The reason is placement.  That translate sits two units from h, so the critical
structure constrains its own neighbourhood and reaches the hexagon only through
the loose ball between them.  What rung one needs is the hexagon's six vertices
to be inside the strain, not beside it:

    no 5-colouring makes both inscribed sqrt3-triangles monochromatic

is a statement about those six points, so each of them should be a vertex of a
5-critical graph in its own right.

So translate the 803 six times, once per hexagon vertex, placing its densest
vertex exactly there.  The copies overlap each other -- the hexagon's vertices
are a unit apart and the copies are far wider than that -- so this is one graph,
not six, and every triangle vertex carries a critical graph's worth of
constraint.  The unit circle about h is pruned to the hexagon alone, so the
burden stays at twenty.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(P)
adj0 = defaultdict(set)
for x, y in g0.edges(): adj0[x].add(y); adj0[y].add(x)
v0 = max(range(g0.n), key=lambda v: len(adj0[v]))
r3 = F.sqrt(3); half = F.rational(Fr(1, 2))
c60, s60 = half, r3 * half
O = Point(F.rational(0), F.rational(0))
hexa = []
x, y = F.rational(1), F.rational(0)
for _ in range(6):
    hexa.append(Point(x, y))
    x, y = x * c60 - y * s60, x * s60 + y * c60
pts = {(0.0, 0.0): O}
for hv in hexa:
    pts[(round(hv.fx, 9), round(hv.fy, 9))] = hv
for hv in hexa:
    tx, ty = hv.x - g0.vertices[v0].x, hv.y - g0.vertices[v0].y
    for p in P:
        r = Point(p.x + tx, p.y + ty)
        pts[(round(r.fx, 9), round(r.fy, 9))] = r
one = F.rational(1)
hexkeys = {(round(hv.fx, 9), round(hv.fy, 9)) for hv in hexa}
keep = {}
for k, p in pts.items():
    dd = p.fx * p.fx + p.fy * p.fy
    if abs(dd - 1.0) < 1e-9 and k not in hexkeys: continue
    keep[k] = p
G = build_graph(list(keep.values())); n = G.n
E = list(G.edges())
adj = defaultdict(set)
for x, y in E: adj[x].add(y); adj[y].add(x)
h = next(i for i in range(n) if (G.vertices[i] - O).norm2() == F.rational(0))
nb = sorted(adj[h])
degs = sorted(len(adj[v]) for v in nb)
print(f"six critical copies on the hexagon: n={n} edges={len(E)} degree "
      f"{2*len(E)/n:.2f}; |N(h)|={len(nb)}, hexagon degrees {degs}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def colourable(k, budget=None):
    X = lambda v, c: 1 + v * k + c
    cnf = [[X(v, c) for c in range(k)] for v in range(n)]
    for x, y in E:
        for c in range(k):
            cnf.append([-X(x, c), -X(y, c)])
    tri = None
    for u in range(n):
        for v in sorted(adj[u]):
            w = adj[u] & adj[v]
            if w: tri = (u, v, min(w)); break
        if tri: break
    if tri:
        for i, v in enumerate(tri): cnf.append([X(v, i)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if budget is None: r = s.solve()
    else:
        s.conf_budget(budget); r = s.solve_limited()
    s.delete(); return r

print(f"  4-colourable: {colourable(4, 30_000_000)}   [{time.time()-t0:.0f}s]",
      flush=True)
r5 = colourable(5, 30_000_000)
print(f"  5-colourable: {r5}   [{time.time()-t0:.0f}s]", flush=True)
if r5 is False:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"source": NAME, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/hexcrit_hit.json", "w"))
    sys.exit(0)
ms = MuSolver(G, 5, budget=None)
r = ms.at_most_two(nb)
print(f"  at_most_two(hexagon) = {r}   "
      f"({'placeable' if r else '*** RUNG ONE ***'})   [{time.time()-t0:.0f}s]",
      flush=True)
if r is False:
    print(f"  *** mu_5 = {ms.mu(nb)} ***", flush=True)
    json.dump({"source": NAME, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/hexcrit_rung1.json", "w"))
ms.close()
