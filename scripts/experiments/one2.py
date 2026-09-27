"""Does this graph map to K(p/q)?  With the symmetries taken out.

Two changes over the naive encoding, both of which matter for UNSAT, which
is the only outcome that proves anything:

At most one position per vertex is not needed.  If an assignment gives a
vertex several positions and every edge clause still holds, then picking any
one position per vertex gives a genuine homomorphism -- the constraints are
all negative, so a smaller assignment can only satisfy more of them.  That
drops p(p-1)/2 clauses per vertex.

K(p/q) has an automorphism group containing all p rotations and a
reflection, so every homomorphism comes in a family of 2p, and a solver
proving UNSAT reproves it 2p times over.  Pinning one vertex to position 0
kills the rotations; confining one of its neighbours to the half of its
allowed range at or below p/2 kills the reflection, which fixes 0 and maps
that range onto itself.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, gc
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
CAR, R = sys.argv[1], Fr(*map(int, sys.argv[2].split("/")))
t0 = time.time()
P = (build_G(K, as_graph=False) if CAR == "G"
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
assert b.overflow_headroom(rows) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
del rows
gc.collect()
n, p, q = len(P), R.numerator, R.denominator
deg = defaultdict(int)
nbr = defaultdict(list)
for a, c in E:
    deg[a] += 1
    deg[c] += 1
    nbr[a].append(c)
    nbr[c].append(a)
v0 = max(range(n), key=lambda v: deg[v])
u0 = max(nbr[v0], key=lambda u: deg[u])
print(f"{CAR}: {n} pts, {len(E)} edges -> K({p}/{q}) = {float(R):.4f}; "
      f"pin v{v0} (deg {deg[v0]}) at 0, confine v{u0} to [{q},{p//2}]"
      f"  [{time.time()-t0:.0f}s]", flush=True)
s = Solver(name="cd15")
for v in range(n):
    s.add_clause([1 + v * p + j for j in range(p)])
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            s.add_clause([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
for j in range(1, p):
    s.add_clause([-(1 + v0 * p + j)])
s.add_clause([1 + v0 * p])
s.add_clause([1 + u0 * p + j for j in range(q, p // 2 + 1)])
print(f"clauses in, solving  [{time.time()-t0:.0f}s]", flush=True)
ok = s.solve()
print(("MAPS  -> chi_c <= " + str(R)) if ok else
      ("DOES NOT MAP  -> chi_c > " + str(R) + f" = {float(R):.4f}"),
      f" [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
