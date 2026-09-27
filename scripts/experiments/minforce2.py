"""Minimise the forcing inside a saved union, with randomised passes and kicks.

The 5-chromatic graph is twice the size of whatever still forces the pair, and
the forcing instance -- "colour this subgraph with four colours, with c(a)=0
and c(b)=1" -- is heavily constrained, so it proves UNSAT in a second or two
where plain 4-colourability of the doubled graph takes minutes.  Minimising
the ingredient is the cheap direction by two orders of magnitude.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, rotation_joining
from hn.graph import build_graph
from pysat.solvers import Solver

SRC, SEED, ROUNDS, KICK = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
d = json.load(open(SRC))
F = Field(tuple(d["field_generators"]))
U = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
A_, B_ = d["forced_pairs"][0]
g = build_graph(U)
n, K = g.n, 4
d2 = (U[A_] - U[B_]).norm2()
print(f"union {n} points; forcing pair ({A_},{B_}) at d^2={d2}", flush=True)
X = lambda v, c: 1 + v * K + c
SEL = lambda v: 1 + n * K + v
cnf = [[-SEL(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
cnf.append([X(A_, 0)]); cnf.append([X(B_, 1)])

def forces(sub):
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[SEL(v) for v in sorted(sub)])
    if ok:
        s.delete(); return None
    core = s.get_core(); s.delete()
    return (set(l - 1 - n * K for l in core) | {A_, B_}) if core else set(sub)

rot = rotation_joining(Fr(d2.c[0]), F).about(U[A_])
def spindle_size(sub):
    seen, Z = set(), []
    for v in sorted(sub):
        for z in (U[v], rot(U[v])):
            if z not in seen: seen.add(z); Z.append(z)
    return len(Z), Z

rng = random.Random(SEED)
cur = set(range(n)); best = set(cur); dropped = set()
t0 = time.time()
assert forces(cur) is not None
for rd in range(ROUNDS):
    c = forces(cur)
    if c is not None and len(c) < len(cur):
        dropped |= (cur - c); cur = c
    order = list(cur); rng.shuffle(order)
    done = 0
    for v in order:
        if v in (A_, B_) or v not in cur: continue
        trial = cur - {v}
        c = forces(trial)
        done += 1
        if c is not None:
            dropped |= (cur - c) if len(c) < len(trial) else {v}
            cur = c if len(c) < len(trial) else trial
    sz, _ = spindle_size(cur)
    print(f"  seed {SEED} round {rd}: forcing {len(cur)} -> spindle {sz}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if len(cur) < len(best):
        best = set(cur)
        sz, Z = spindle_size(best)
        json.dump({"field_generators": list(F.gens), "n_forcing": len(best),
                   "n": sz, "seed": SEED, "source": SRC,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in Z]},
                  open(f"/tmp/hn/mf2_{SEED}.json", "w"))
    if dropped:
        cur = cur | set(rng.sample(sorted(dropped), min(KICK, len(dropped))))
sz, _ = spindle_size(best)
print(f"seed {SEED}: best forcing {len(best)}, spindle {sz}   [{time.time()-t0:.0f}s]",
      flush=True)
