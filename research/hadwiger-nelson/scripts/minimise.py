"""Minimise a 5-chromatic unit-distance graph, with kicks.

Pure greedy deletion stops at the first subset that is minimal under single
deletions, which is usually far from small.  The standard escape is to perturb:
once greedy converges, put back a random handful of deleted vertices and run
greedy again from a different order.  Each round is cheap because the instance
gets smaller, and a round that does not improve costs nothing but time.

Accepts a starting point set from JSON (the exact field-basis vectors) so it
can be pointed at Z, at Z' from the forcing shrink, or at any later carrier.
"""
import sys, time, json, random
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

SRC = sys.argv[1]
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
KICK = int(sys.argv[3]) if len(sys.argv) > 3 else 25
OUT = sys.argv[4] if len(sys.argv) > 4 else "/home/user/darwin-50/research/hadwiger-nelson/data/five_247_min.json"

d = json.load(open(SRC))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
K = 4
print(f"start: {len(P)} points in Q{F.gens}", flush=True)

def refuses(sub):
    """Does the induced subgraph on `sub` refuse four colours?  Returns an
    UNSAT core (a smaller refusing subset) or None."""
    g = build_graph([P[i] for i in sub])
    n = g.n
    X = lambda v, c: 1 + v * K + c
    A = lambda v: 1 + n * K + v
    cnf = [[-A(v)] + [X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[A(v) for v in range(n)])
    if ok:
        s.delete(); return None
    core = s.get_core(); s.delete()
    if not core: return list(sub)
    keep = set(l - 1 - n * K for l in core)
    return [sub[i] for i in sorted(keep)]

t0 = time.time()
cur = list(range(len(P)))
assert refuses(cur) is not None, "the input is 4-colourable"
best = cur[:]
rng = random.Random(2026)
dropped = []
for rd in range(ROUNDS):
    c = refuses(cur)
    if c is not None and len(c) < len(cur):
        dropped += [v for v in cur if v not in set(c)]
        cur = c
    order = cur[:]
    rng.shuffle(order)
    for v in order:
        if v not in cur: continue
        trial = [u for u in cur if u != v]
        c = refuses(trial)
        if c is not None:
            dropped += [u for u in cur if u not in set(c)]
            cur = c
    print(f"  round {rd}: {len(cur)} vertices   [{time.time()-t0:.0f}s]", flush=True)
    if len(cur) < len(best):
        best = cur[:]
        g = build_graph([P[i] for i in best])
        json.dump({"field_generators": list(F.gens), "n": g.n, "source": SRC,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in g.vertices]}, open(OUT, "w"))
        print(f"    -> new best, written to {OUT}", flush=True)
    # kick: put a handful of deleted vertices back and let greedy re-decide
    if dropped:
        add = rng.sample(dropped, min(KICK, len(dropped)))
        cur = sorted(set(cur) | set(add))
g = build_graph([P[i] for i in best])
print(f"\nbest: n={g.n} m={sum(len(a) for a in g.adj)//2}   [{time.time()-t0:.0f}s]",
      flush=True)
