"""Minimise the forcing subgraph, with randomised passes and kicks.

Z refuses four colours because H forces one pair to share a colour and the
rotation puts the two images of the far endpoint one apart.  So the size of Z
is governed by the size of the smallest subgraph that still forces the pair,
and THAT instance -- "colour this subgraph with four colours, with c(a)=0 and
c(b)=1" -- is heavily constrained and proves UNSAT in seconds, where "colour
Z with four colours" takes minutes.  Minimising the ingredient is the cheap
direction by two orders of magnitude.

One greedy pass stops at the first subset minimal under single deletions,
which is rarely small.  Each round here shuffles the deletion order, and when
a round stops improving a handful of deleted vertices go back so the next
round starts somewhere else.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

SEED = int(sys.argv[1]); ROUNDS = int(sys.argv[2]); KICK = int(sys.argv[3])
K1 = Field((3, 11, 247))
pts = build_Sa(K1)
r1 = rotation_joining(Fr(1), K1).about(pts[25])
seen, H = set(), []
for p in pts:
    for q in (p, r1(p)):
        if q not in seen: seen.add(q); H.append(q)
g = build_graph(H)
n, K = g.n, 4
A_, B_ = 157, 327
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

rng = random.Random(SEED)
cur = set(range(n)); best = set(cur); dropped = set()
t0 = time.time()
for rd in range(ROUNDS):
    c = forces(cur)
    if c is not None and len(c) < len(cur):
        dropped |= (cur - c); cur = c
    order = list(cur); rng.shuffle(order)
    for v in order:
        if v in (A_, B_) or v not in cur: continue
        trial = cur - {v}
        c = forces(trial)
        if c is not None:
            dropped |= (cur - c) if len(c) < len(trial) else {v}
            cur = c if len(c) < len(trial) else trial
    print(f"  seed {SEED} round {rd}: {len(cur)}   [{time.time()-t0:.0f}s]", flush=True)
    if len(cur) < len(best):
        best = set(cur)
        rot = rotation_joining(Fr(64, 9), K1).about(H[A_])
        Hp = [H[v] for v in sorted(best)]
        seen2, Zp = set(), []
        for p in Hp:
            for q in (p, rot(p)):
                if q not in seen2: seen2.add(q); Zp.append(q)
        gz = build_graph(Zp)
        json.dump({"field_generators": list(K1.gens), "n_forcing": len(best),
                   "n": gz.n, "seed": SEED,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gz.vertices]},
                  open(f"/tmp/hn/minf_{SEED}.json", "w"))
        print(f"    -> H'={len(best)}, Z'={gz.n} points saved", flush=True)
    if dropped:
        cur = cur | set(rng.sample(sorted(dropped), min(KICK, len(dropped))))
print(f"seed {SEED} best H' = {len(best)}   [{time.time()-t0:.0f}s]", flush=True)
