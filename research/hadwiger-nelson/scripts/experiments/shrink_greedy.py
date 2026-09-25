"""Hunt a small 5-chromatic graph: many greedy orders, and a budget per test.

Order is everything.  "Is G - v still 4-uncolourable?" is FAST when the answer
is no -- v is essential, G - v is 4-colourable, and a colouring comes straight
back -- and SLOW when the answer is yes, since a removal has to be certified by
a refutation.  So a pass costs about one refutation per vertex actually removed,
and greedy stops at whatever vertex-critical subgraph its order walks into.  The
803 came from one such walk; a different order lands somewhere else, which is
why this is a search and not a computation.

Two levers.  Random restarts, keeping the best.  And starting from the DENSE
carriers rather than the sparse ones: the same refusal of four carried by more
edges may need far fewer vertices, and this project's dense 5-chromatic graphs
run to mean degree 18.5 where the 803 has 10.1.

A conflict budget keeps one hard refutation from eating the pass: a test that
runs out is treated as "essential", which can only make the answer larger, never
wrong.
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
SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
BUDGET = int(sys.argv[3]) if len(sys.argv) > 3 else 2_000_000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
print(f"{NAME} n={n} edges={len(E)} mean degree {2*len(E)/n:.2f}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c

def refuses_four(keep, budget=BUDGET):
    """True if G[keep] has no 4-colouring.  None if the budget ran out."""
    ks = set(keep)
    cnf = [[X(v, c) for c in range(K)] for v in keep]
    for x, y in E:
        if x in ks and y in ks:
            for c in range(K):
                cnf.append([-X(x, c), -X(y, c)])
    # pin a triangle inside `keep` to kill the 4! symmetric copies of every
    # refutation; sound for deciding colourability, and worth a factor of sixty
    tri = None
    for u in keep:
        for v in sorted(adj[u] & ks):
            w = adj[u] & adj[v] & ks
            if w:
                tri = (u, v, min(w)); break
        if tri: break
    if tri:
        for i, v in enumerate(tri):
            cnf.append([X(v, i)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if budget is None:
        r = s.solve()
    else:
        s.conf_budget(budget); r = s.solve_limited()
    s.delete()
    if r is None: return None
    return not r

base = list(range(n))
if refuses_four(base, None) is not True:
    print("  this graph does not refuse four", flush=True); sys.exit(0)
print(f"  refuses four   [{time.time()-t0:.0f}s]", flush=True)
best = None
for seed in range(SEEDS):
    cur = list(base)
    order = list(base); random.seed(seed); random.shuffle(order)
    removed = 0
    for j, v in enumerate(order):
        if v not in cur: continue
        trial = [w for w in cur if w != v]
        r = refuses_four(trial)
        if r is True:
            cur = trial; removed += 1
            if removed % 50 == 0:
                print(f"    seed {seed}: {len(cur)} left ({removed} removed)   "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    m = sum(1 for x, y in E if x in set(cur) and y in set(cur))
    print(f"  seed {seed}: critical subgraph {len(cur)} vertices, {m} edges   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if best is None or len(cur) < best[0]:
        best = (len(cur), cur)
        json.dump({"source": NAME, "seed": seed, "n": len(cur),
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in g.vertices[v].x.c],
                               [[t.numerator, t.denominator] for t in g.vertices[v].y.c]]
                              for v in cur]},
                  open(f"{ROOT}/data/hunt_{len(cur)}.json", "w"))
        print(f"    written data/hunt_{len(cur)}.json", flush=True)
print(f"\n  best: {best[0]} vertices   [{time.time()-t0:.0f}s]", flush=True)
