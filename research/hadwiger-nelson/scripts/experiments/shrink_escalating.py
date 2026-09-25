"""Escalating budgets: take the cheap removals first, they pay for the rest.

The removal test is fast when the vertex is essential -- a colouring comes
straight back -- and slow when it is removable, because a removal has to be
certified by a refutation.  So a single pass at a generous budget spends all its
time on the removals, forty-seven seconds each on the 1139-point graph, and
twenty-five of them is twenty minutes before anything is even reported.

Escalate instead.  Sweep at a small conflict budget first: the easy removals
land immediately, each one shrinking the instance and making every later call
cheaper, and a test that runs out of budget is simply treated as essential --
which can only leave the answer larger, never wrong.  Then sweep again with more
budget on what survived, and again.  Later sweeps succeed where earlier ones
timed out, because the graph they are asked about is smaller.

Order by degree ascending, for the same reason: a low-degree vertex constrains
least, so it is the likeliest to be removable and the cheapest to certify.

One warm solver throughout, with selectors, so nothing is ever rebuilt, and the
UNSAT core after each success drops every vertex it did not need.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 4
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247.json"
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 0
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
tri = None
for u in range(n):
    for v in sorted(adj[u]):
        w = adj[u] & adj[v]
        if w: tri = (u, v, min(w)); break
    if tri: break
print(f"{NAME} n={n} edges={len(E)} degree {2*len(E)/n:.2f}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
S = lambda v: n * K + 1 + v
cnf = [[-S(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])
for i, v in enumerate(tri):
    cnf.append([-S(v), X(v, i)])
solver = Solver(name="cd19", bootstrap_with=cnf)
PINNED = {S(v) for v in tri}

def refuses(alive, budget):
    solver.conf_budget(budget)
    r = solver.solve_limited(assumptions=sorted(alive))
    if r is None or r: return False, None
    return True, set(solver.get_core() or alive)

ok, core = refuses({S(v) for v in range(n)}, 50_000_000)
if not ok:
    print("  does not refuse four", flush=True); sys.exit(0)
alive = set(core) | PINNED
print(f"  refuses four; first core {len(alive)}   [{time.time()-t0:.0f}s]",
      flush=True)
random.seed(SEED)
best = len(alive)
# The small budgets were wasted.  A removal is certified by a REFUTATION, and
# twenty thousand conflicts never reaches one on a near-critical graph, so the
# cheap sweep spends forty minutes proving nothing while paying two seconds a
# time for the essential vertices it cannot skip.  Start where a removal can
# actually land.
for budget in (2_000_000, 8_000_000, 30_000_000):
    changed = True
    while changed:
        changed = False
        order = sorted((lit for lit in alive if lit not in PINNED),
                       key=lambda lit: (len(adj[lit - (n * K + 1)]),
                                        random.random()))
        tested = 0
        for lit in order:
            if lit not in alive: continue
            ok, c2 = refuses(alive - {lit}, budget)
            tested += 1
            if ok:
                alive = (c2 & (alive - {lit})) | PINNED
                changed = True
            if tested % 100 == 0:
                print(f"      budget {budget}: {tested} tested, {len(alive)} "
                      f"alive   [{time.time()-t0:.0f}s]", flush=True)
        if len(alive) < best:
            best = len(alive)
            verts = sorted(a - (n * K + 1) for a in alive)
            m = sum(1 for x, y in E if x in set(verts) and y in set(verts))
            print(f"    budget {budget}: {len(verts)} vertices, {m} edges, "
                  f"degree {2*m/len(verts):.2f}   [{time.time()-t0:.0f}s]",
                  flush=True)
            json.dump({"source": NAME, "seed": SEED, "n": len(verts),
                       "field_generators": list(F.gens),
                       "points": [[[[t.numerator, t.denominator]
                                    for t in g.vertices[v].x.c],
                                   [[t.numerator, t.denominator]
                                    for t in g.vertices[v].y.c]]
                                  for v in verts]},
                      open(f"{ROOT}/data/hunt_{len(verts)}.json", "w"))
print(f"\n  final: {best} vertices   [{time.time()-t0:.0f}s]", flush=True)
