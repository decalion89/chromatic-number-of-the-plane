"""Frozen literals: a progress measure that looks at the space, not the sample.

The kill rate was rejected for good reason -- it counts how easy the sampled
homomorphisms are to block, and the sample is a vanishing part of an
astronomical set.  Rejecting it left the hole unfilled, and this fills it.

For a vertex v and a position j, ask the solver for a homomorphism with v at
j.  UNSAT means no homomorphism anywhere puts v there: the literal is frozen.
Count the frozen pairs and you have measured the solution space itself, not a
draw from it.  The count only rises as points are added, and it reaches n*p
exactly when the graph refuses the ratio outright, so it is a real distance
to the goal rather than a proxy for one.

Assumption-based incremental solving makes it affordable: one solver, one
assumption per query, all the learnt clauses kept between queries.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR, R = sys.argv[1], Fr(*map(int, sys.argv[2].split("/")))
SAMPLE = int(sys.argv[3]) if len(sys.argv) > 3 else 60
t0 = time.time()
P = (build_G(K, as_graph=False) if CAR == "G"
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), R.numerator, R.denominator
deg = defaultdict(int)
for a, c in E:
    deg[a] += 1
    deg[c] += 1
s = Solver(name="cd15")
for v in range(n):
    s.add_clause([1 + v * p + j for j in range(p)])
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            s.add_clause([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
print(f"{CAR}: {n} pts, {len(E)} edges, ratio {R} = {float(R):.4f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
if not s.solve():
    print("graph already REFUSES the ratio -- every literal is frozen",
          flush=True)
    sys.exit()
# pin one vertex to kill the p rotations, else nothing is ever frozen
v0 = max(range(n), key=lambda v: deg[v])
for j in range(1, p):
    s.add_clause([-(1 + v0 * p + j)])
s.add_clause([1 + v0 * p])
random.seed(7)
vs = sorted(random.sample([v for v in range(n) if v != v0],
                          min(SAMPLE, n - 1)))
froz = tot = 0
for k, v in enumerate(vs):
    for j in range(p):
        tot += 1
        if not s.solve(assumptions=[1 + v * p + j]):
            froz += 1
    if k % 15 == 14:
        print(f"   {k+1}/{len(vs)} vertices: {froz}/{tot} frozen "
              f"({100*froz/tot:.1f}%)  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nFROZEN {froz}/{tot} = {100*froz/tot:.2f}% of sampled literals"
      f"   (100% would mean the graph refuses {R})"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
