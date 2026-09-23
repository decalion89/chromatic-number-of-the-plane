"""Forced pairs at five, sought where the pricing theorem says they must live.

Two results from this session, put together.

The pricing: a forced pair's certificate is the WHOLE near-critical graph --
every genuine forced pair at four colours needed all 802 vertices of the
vertex-critical 803 minus one.  So forced pairs are not local objects; they
exist only in graphs sitting at the edge of the next chromatic number, and
hunting them in loose 5-chromatic graphs was always going to come back empty.
It did, over twenty-five million pairs.

The tightness: six copies of the vertex-critical 803, each placed with its
densest vertex on one vertex of a hexagon and each TURNED by a different
multiple of 60 degrees, give a 4159-point graph that refuses four and admits
five -- but the solver needed 452 seconds to find the five-colouring, where the
803 takes under one.  Nothing built here has sat that close to the edge.

So look for forced pairs HERE.  And pose the target accordingly: not a graph with
no 5-colouring, but one with essentially ONE.  In a uniquely 5-colourable graph
every non-adjacent pair is forced, together or apart, and a single forced-
together pair at any distance d in (0,2), d != 1, is closed by one spindle.

Method, as in the 25-million-pair settlement: one warm solver, each candidate
pair tested by assuming it split, and every colouring that comes back used to
eliminate every other candidate it splits.  Candidates restricted to the pairs a
spindle can use, 0 < d < 2 and d != 1.
"""
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = "five_247_c.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(P)
adj0 = defaultdict(set)
for x, y in g0.edges(): adj0[x].add(y); adj0[y].add(x)
v0 = max(range(g0.n), key=lambda v: len(adj0[v]))
r3 = F.sqrt(3); half = F.rational(Fr(1, 2)); c60, s60 = half, r3 * half
O = Point(F.rational(0), F.rational(0))
hexa = []; x, y = F.rational(1), F.rational(0)
for _ in range(6):
    hexa.append(Point(x, y)); x, y = x * c60 - y * s60, x * s60 + y * c60
pts = {(0.0, 0.0): O}
for hv in hexa: pts[(round(hv.fx, 9), round(hv.fy, 9))] = hv
for i, hv in enumerate(hexa):
    cc, ss = F.rational(1), F.rational(0)
    for _ in range(i): cc, ss = cc * c60 - ss * s60, cc * s60 + ss * c60
    bx, by = g0.vertices[v0].x, g0.vertices[v0].y
    for p in P:
        dx, dy = p.x - bx, p.y - by
        r = Point(hv.x + dx * cc - dy * ss, hv.y + dx * ss + dy * cc)
        pts[(round(r.fx, 9), round(r.fy, 9))] = r
hexkeys = {(round(hv.fx, 9), round(hv.fy, 9)) for hv in hexa}
keep = {k: p for k, p in pts.items()
        if not (abs(p.fx * p.fx + p.fy * p.fy - 1.0) < 1e-9 and k not in hexkeys)}
G = build_graph(list(keep.values())); n = G.n
E = list(G.edges())
adj = defaultdict(set)
for a, b in E: adj[a].add(b); adj[b].add(a)
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
print(f"tight graph: n={n} edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"construction": "six turned copies of five_247_c on a hexagon",
           "field_generators": list(F.gens),
           "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                       [[t.numerator, t.denominator] for t in q.y.c]]
                      for q in G.vertices]},
          open(f"{ROOT}/data/tight_hexagon_4159.json", "w"))

X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            cnf.append([-X(v, a), -X(v, b)])
for a, b in E:
    for c in range(K):
        cnf.append([-X(a, c), -X(b, c)])
s = Solver(name="cd19", bootstrap_with=cnf)
if not s.solve():
    print("  *** NOT 5-COLOURABLE ***", flush=True); sys.exit(0)
def colouring():
    m = s.get_model()
    return [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
col = colouring()
print(f"  first colouring   [{time.time()-t0:.0f}s]", flush=True)
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
cand = []
for i in range(n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i or j in adj[i]: continue
                dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
                if 1e-12 < dd < 4.0 - 1e-9 and col[i] == col[j]:
                    cand.append((i, j))
print(f"  {len(cand)} spindle-able pairs monochromatic in the first colouring   "
      f"[{time.time()-t0:.0f}s]", flush=True)
nsel = n * K + 1
forced = []
rounds = 0
while cand:
    i, j = cand[0]
    sel = nsel; nsel += 1
    for c in range(K):
        s.add_clause([-sel, -X(i, c), -X(j, c)])
    r = s.solve(assumptions=[sel])
    s.add_clause([-sel])                      # retire the selector
    if not r:
        forced.append((i, j))
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        print(f"  *** FORCED TOGETHER AT FIVE: ({i},{j}), d^2 = {dd:.6f} ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        json.dump({"graph": "tight_hexagon_4159.json", "forced": forced},
                  open(f"{ROOT}/data/tight_forced.json", "w"))
        cand = cand[1:]
        continue
    c2 = colouring(); rounds += 1
    cand = [(a, b) for a, b in cand if c2[a] == c2[b]]
    if rounds % 5 == 0 or len(cand) < 50:
        print(f"    colouring {rounds}: {len(cand)} candidates left   "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\n  {len(forced)} forced pairs; {rounds} colourings used   "
      f"[{time.time()-t0:.0f}s]", flush=True)
