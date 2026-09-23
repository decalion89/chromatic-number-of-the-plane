"""Six critical copies on the hexagon, each turned a different way.

The translated version -- six copies of the vertex-critical 803, one per hexagon
vertex, all in the same orientation -- gives 3830 points, degree 12.69, hexagon
degrees all 30, and refuses four.  But copies related by pure translation are
"parallel": whatever structure one imposes, its neighbour imposes the same thing
shifted, and the union can inherit a colouring that respects the shift.

Turning each copy breaks that.  Rotate copy i about its own hexagon vertex --
which keeps that vertex exactly where it must be -- by i times a field angle, so
the six critical structures interlock at six different attitudes and no
translation of a colouring carries between them.  60 degrees needs only sqrt3
and is already in the field; the spindle angle of the carrier is another choice
and is less symmetric, which may matter more.

Everything else is unchanged and is what the theory asked for: N(h) pruned to a
single hexagon, so one component and twenty patterns, and every hexagon vertex
the densest vertex of a 5-critical graph in its own right.
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
MODE = sys.argv[2] if len(sys.argv) > 2 else "turn"
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
    hexa.append(Point(x, y)); x, y = x * c60 - y * s60, x * s60 + y * c60
pts = {(0.0, 0.0): O}
for hv in hexa: pts[(round(hv.fx, 9), round(hv.fy, 9))] = hv
for i, hv in enumerate(hexa):
    cc, ss = F.rational(1), F.rational(0)
    for _ in range(i if MODE == "turn" else 0):
        cc, ss = cc * c60 - ss * s60, cc * s60 + ss * c60
    bx, by = g0.vertices[v0].x, g0.vertices[v0].y
    for p in P:
        dx, dy = p.x - bx, p.y - by
        rx, ry = dx * cc - dy * ss, dx * ss + dy * cc
        r = Point(hv.x + rx, hv.y + ry)
        pts[(round(r.fx, 9), round(r.fy, 9))] = r
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
print(f"mode={MODE}: n={n} edges={len(E)} degree {2*len(E)/n:.2f}; "
      f"|N(h)|={len(nb)}, hexagon degrees {sorted(len(adj[v]) for v in nb)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def colourable(k, budget):
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
    s.conf_budget(budget); r = s.solve_limited(); s.delete(); return r

print(f"  4-colourable: {colourable(4, 30_000_000)}   [{time.time()-t0:.0f}s]",
      flush=True)
r5 = colourable(5, 30_000_000)
print(f"  5-colourable: {r5}   [{time.time()-t0:.0f}s]", flush=True)
if r5 is False:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"source": NAME, "mode": MODE, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/hexcrit2_hit.json", "w"))
    sys.exit(0)
ms = MuSolver(G, 5, budget=None)
r = ms.at_most_two(nb)
print(f"  at_most_two(hexagon) = {r}   "
      f"({'placeable' if r else '*** RUNG ONE ***'})   [{time.time()-t0:.0f}s]",
      flush=True)
if r is False:
    print(f"  *** mu_5 = {ms.mu(nb)} ***", flush=True)
    json.dump({"source": NAME, "mode": MODE, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/hexcrit2_rung1.json", "w"))
ms.close()
