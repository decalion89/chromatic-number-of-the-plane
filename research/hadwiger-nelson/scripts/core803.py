"""Minimise the 803-point 5-chromatic graph by UNSAT cores, not by vertices.

Everything at five colours is priced by the size of the carrier: one
5-colouring of the 7141-point symmetric graph costs cadical 322 s and of the
13261-point one 1718 s, so a filter needing hundreds of them is affordable
only on a small object.  The record for a 5-chromatic unit-distance graph is
the 509-vertex Parts graph, reached from de Grey's 1581 by minimisations
costing on the order of 100 000 CPU-hours.  This one is 803, in a different
field, and has only ever been cut by removing whole orbits.

Greedy vertex removal needs one 4-colour UNSAT proof per vertex, which is days.
An UNSAT core gives a large cut per solve instead.  Gate every constraint of
vertex v behind a selector s_v -- the at-least-one clause becomes
[-s_v, X(v,0..K-1)] and each edge clause [-s_x, -s_y, ...] -- so switching s_v
off removes v and all its edges, exactly the induced subgraph.  Assume every
selector true; the solver returns UNSAT with a core, and that core is a set of
vertices that already refuses four.  Restrict to it and repeat until the core
stops shrinking, then finish greedily on what is left.

The graph stays a unit-distance graph throughout -- vertices are only deleted,
never moved -- so the field and the certificate survive untouched.
"""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 247))
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
PTS = [Point(F.element([Fr(a, b) for a, b in x]),
             F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
K = 4          # refusing FOUR is what makes it 5-chromatic

def edges_of(pts):
    g = build_graph(pts)
    return g, list(g.edges())

def core_round(keep):
    """One core extraction on the induced subgraph over `keep`."""
    pts = [PTS[i] for i in keep]
    g, E = edges_of(pts)
    n = g.n
    X = lambda v, c: 1 + v * K + c
    sel = lambda v: 1 + n * K + v
    cl = [[-sel(v)] + [X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in E:
        for c in range(K):
            cl.append([-sel(x), -sel(y), -X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cl)
    r = s.solve(assumptions=[sel(v) for v in range(n)])
    if r:
        s.delete()
        return None, n, len(E)
    core = s.get_core(); s.delete()
    back = {sel(v): v for v in range(n)}
    idx = sorted(back[l] for l in core if l in back)
    return [keep[i] for i in idx], n, len(E)

keep = list(range(len(PTS)))
g0, E0 = edges_of(PTS)
print(f"start n={g0.n} m={len(E0)}", flush=True)
while True:
    nxt, n, m = core_round(keep)
    if nxt is None:
        print(f"  !! the subgraph became 4-colourable at n={n} -- stop", flush=True)
        break
    print(f"  core: {n} -> {len(nxt)} vertices   [{time.time()-t0:.0f}s]", flush=True)
    if len(nxt) >= len(keep):
        keep = nxt
        break
    keep = nxt

# greedy finish: drop a vertex whenever the rest still refuses four
pts = [PTS[i] for i in keep]
g, E = edges_of(pts)
print(f"  after cores: n={g.n} m={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

def refuses_four(sub):
    q = [PTS[i] for i in sub]
    gg = build_graph(q); nn = gg.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(nn)]
    for v in range(nn):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in gg.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    sv = Solver(name="cd19", bootstrap_with=cl); r = sv.solve(); sv.delete()
    return not r

order = sorted(range(len(keep)),
               key=lambda i: len(build_graph([PTS[j] for j in keep])
                                 .adj[i]))
cur = list(keep)
removed = 0
for step, v in enumerate(list(cur)):
    trial = [u for u in cur if u != v]
    if len(trial) < 4:
        continue
    if refuses_four(trial):
        cur = trial; removed += 1
    if step % 25 == 0:
        print(f"    greedy {step}/{len(keep)}: n={len(cur)} "
              f"({removed} dropped)   [{time.time()-t0:.0f}s]", flush=True)
        gg = build_graph([PTS[i] for i in cur])
        json.dump({"field_generators": list(F.gens), "n": gg.n,
                   "m": sum(len(a) for a in gg.adj) // 2,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gg.vertices]},
                  open(f"{ROOT}/data/five_247_min.json", "w"))
gg = build_graph([PTS[i] for i in cur])
mm = sum(len(a) for a in gg.adj) // 2
print(f"  FINAL n={gg.n} m={mm} deg={2.0*mm/gg.n:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
json.dump({"field_generators": list(F.gens), "n": gg.n, "m": mm,
           "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                       [[c.numerator, c.denominator] for c in p.y.c]]
                      for p in gg.vertices]},
          open(f"{ROOT}/data/five_247_min.json", "w"))
print("  written data/five_247_min.json", flush=True)
