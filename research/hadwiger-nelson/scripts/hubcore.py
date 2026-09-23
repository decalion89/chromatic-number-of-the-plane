"""The narrowest hub disjunction, found by a core rather than by rings.

The consumable shape is a statement about ONE point:

    every 5-colouring gives  c(h) = c(v)  for some v in S,

because rotations about h fix h, so every rotated copy says this about the same
h, and all the partners it can name are the rotated images of S -- which stay on
the same circles about h, where the clash geometry is computable.  Nothing
requires S to be a whole ring; what matters is that S be SMALL.

The statement holds at all (with S = everything) exactly when G - h is not
4-colourable, so hubs live outside the graph's 5-critical core -- and this
project already owns one: five_247_c is a vertex-critical 5-chromatic subgraph
of five_247, so every vertex of the 1139 that is not one of its 803 is a hub, no
solving needed.

Then a selector per partner and one UNSAT core does what thousands of solves
would, and a greedy pass makes S minimal.

Pigeonhole-free with no check: every extra edge is incident to h, and a
unit-distance graph in the plane has clique number 3, so no K6 can appear.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
BIG, CORE = "five_247.json", "five_247_c.json"
db = json.load(open(f"{ROOT}/data/{BIG}"))
dc = json.load(open(f"{ROOT}/data/{CORE}"))
F = Field(tuple(db["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in db["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
Fc = Field(tuple(dc["field_generators"]))
Pc = [Point(Fc.element([Fr(a, b) for a, b in x]),
            Fc.element([Fr(a, b) for a, b in y])) for x, y in dc["points"]]
ck = {(round(q.fx, 9), round(q.fy, 9)) for q in Pc}
inside = [i for i in range(n) if (round(hx[i], 9), round(hy[i], 9)) in ck]
outside = [i for i in range(n) if i not in set(inside)]
print(f"{BIG} n={n}; {len(inside)} of the {len(ck)} core points found inside, "
      f"{len(outside)} vertices outside the core   [{time.time()-t0:.0f}s]",
      flush=True)
if len(inside) < len(ck):
    print("  the core is NOT a subgraph -- falling back to solving the filter",
          flush=True)

X = lambda v, c: 1 + v * K + c
S = lambda v: n * K + 1 + v
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])

widths = []
cands = outside if len(inside) == len(ck) else list(range(n))
cands = sorted(cands, key=lambda v: -len(adj[v]))[:15]
for h in cands:
    part = [v for v in range(n) if v != h and v not in adj[h]]
    cnf = list(base)
    for v in part:
        for c in range(K):
            cnf.append([-S(v), -X(h, c), -X(v, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    ass = [S(v) for v in part]
    if s.solve(assumptions=ass):
        print(f"  h={h}: not a hub (G-h is 4-colourable)   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        s.delete(); continue
    core = set(s.get_core() or [])
    for _ in range(40):
        a2 = sorted(core)
        if s.solve(assumptions=a2): break
        c2 = set(s.get_core() or [])
        if len(c2) >= len(core): break
        core = c2
    cur = sorted(core); i = 0
    while i < len(cur):
        trial = cur[:i] + cur[i+1:]
        if trial and not s.solve(assumptions=trial):
            c2 = set(s.get_core() or trial)
            cur = [a for a in trial if a in c2] or trial
            i = 0
        else:
            i += 1
    ps = [a - (n * K + 1) for a in cur]
    ds = Counter(round((hx[h]-hx[v])**2 + (hy[h]-hy[v])**2, 6) for v in ps)
    print(f"  *** hub h={h} deg={len(adj[h])}: |S| = {len(ps)}   "
          f"radii^2 {dict(sorted(ds.items()))}   [{time.time()-t0:.0f}s]", flush=True)
    widths.append({"hub": h, "S": ps, "width": len(ps),
                   "radii2": {str(k): v for k, v in ds.items()}})
    json.dump({"graph": BIG, "hubs": widths},
              open(f"{ROOT}/data/hubcore.json", "w"))
    s.delete()
print(f"\n  {len(widths)} hubs; narrowest "
      f"{min((w['width'] for w in widths), default=None)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
