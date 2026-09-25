"""Shrink with ONE warm solver and a core jump, instead of a fresh CNF each time.

The removal test is "does G - v still refuse four?", and it is slow exactly when
the answer is yes, because a removal has to be certified by a refutation.
Rebuilding the CNF for every test throws away everything the solver learned:
ninety-seven seconds per test on the 6925-point carrier, so fifty removals is an
hour and a pass is a day.

Selectors fix both halves.  Give every vertex s_v, make "v has a colour"
conditional on it, and the test becomes one assumption call on a solver that is
never rebuilt -- an inactive vertex takes no colour and its edges go vacuous, so
the answer is unchanged.  And when the call comes back UNSAT it comes back with
a CORE: the subset of selectors that caused the refutation.  Everything outside
it can be dropped at once, so a single successful removal often takes dozens of
vertices with it instead of one.

A triangle stays pinned to colours 0, 1, 2 and permanently assumed, which kills
the 4! symmetric copies of every refutation.  A conflict budget keeps one hard
call from eating the pass; a test that runs out is treated as essential, which
can only make the answer larger, never wrong.
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
SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 5
BUDGET = int(sys.argv[3]) if len(sys.argv) > 3 else 4_000_000
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
print(f"{NAME} n={n} edges={len(E)} degree {2*len(E)/n:.2f}; triangle {tri}   "
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

def refuses(alive):
    solver.conf_budget(BUDGET)
    r = solver.solve_limited(assumptions=sorted(alive))
    if r is None: return None, None
    if r: return False, None
    return True, set(solver.get_core() or alive)

ok, core = refuses({S(v) for v in range(n)})
print(f"  refuses four: {ok}; first core {len(core) if core else '-'}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
if ok is not True:
    sys.exit(0)
start = set(core) | {S(v) for v in tri}
best = None
for seed in range(SEEDS):
    alive = set(start)
    order = [S(v) for v in range(n) if S(v) in alive and v not in tri]
    random.seed(seed); random.shuffle(order)
    removed = 0
    for lit in order:
        if lit not in alive: continue
        r, c2 = refuses(alive - {lit})
        if r is True:
            alive = (c2 & (alive - {lit})) | {S(v) for v in tri}
            removed += 1
            if removed % 25 == 0:
                print(f"    seed {seed}: {len(alive)} left   "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    verts = sorted(a - (n * K + 1) for a in alive)
    m = sum(1 for x, y in E if x in set(verts) and y in set(verts))
    print(f"  seed {seed}: {len(verts)} vertices, {m} edges, degree "
          f"{2*m/len(verts):.2f}   [{time.time()-t0:.0f}s]", flush=True)
    if best is None or len(verts) < best:
        best = len(verts)
        json.dump({"source": NAME, "seed": seed, "n": len(verts),
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in g.vertices[v].x.c],
                               [[t.numerator, t.denominator] for t in g.vertices[v].y.c]]
                              for v in verts]},
                  open(f"{ROOT}/data/hunt_{len(verts)}.json", "w"))
        print(f"    written data/hunt_{len(verts)}.json", flush=True)
print(f"\n  best: {best}   [{time.time()-t0:.0f}s]", flush=True)
