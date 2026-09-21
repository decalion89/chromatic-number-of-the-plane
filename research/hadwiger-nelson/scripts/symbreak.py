"""Break the colour symmetry with a clique, then ask again.

The instance -- Sa at five colours with classes 4/9 and 16/9 forbidden, 397
points and 13012 clauses -- has held three CDCL solvers for an hour and a
half, and a calibrated local search reaches one conflict and never zero, in
sixty short runs and in a three-million-move deep run.

The formula has a symmetry nobody has taken out of it.  Nothing pins a colour,
so every assignment comes in 120 copies under permuting the five colours, and
a CDCL proof has to rediscover that 120 times over.  The combined graph has
clique number four, so it contains a K4, and any colouring can be permuted to
give that K4 the colours 0, 1, 2, 3 -- fixing them is sound and cuts the
search space by a factor of 120 outright.

Sound is the operative word: the four vertices are pairwise adjacent in the
COMBINED graph, which includes the forbidden pairs as edges, so they take four
distinct colours in any solution of this formula and the permutation exists.

Three solvers again, on the reduced formula.

And a tactical point: prove the EASIER statement first.  The local search
plateaus at one conflict on the two-class instance and at twenty-one on the
four-class one.  If both are unsatisfiable, the four-class one is
unsatisfiable with room to spare and its proof is far shorter -- a
conflict-driven search has to close a much wider gap on the tight one.  The
four-class statement is weaker ("some pair at 2/3, 4/3, 2 or 4") but it is
still a weak property at five colours, which is what nothing here has yet
produced.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k = 5
SETS = {"two": {Fr(4, 9), Fr(16, 9)},
        "three": {Fr(4, 9), Fr(16, 9), Fr(4)},
        "four": {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}}
name = sys.argv[1]
tag = sys.argv[2] if len(sys.argv) > 2 else "two"
CS = SETS[tag]
t0 = time.time()
P = build_Sa(F)
g = build_graph(P)
E = set((min(a, b), max(a, b)) for a, b in g.edges())
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
extra = []
for i in range(len(P) - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in E and Fr(int(sq[off, 0]), D2) in CS:
            extra.append((i, j))
allp = sorted(E) + extra
n = g.n
A = np.zeros((n, n), dtype=bool)
for a, b in allp:
    A[a, b] = A[b, a] = True

# the largest clique reachable greedily in the COMBINED graph.  The first
# attempt took Nc[0] without excluding the vertices already chosen, so it
# could return a "clique" containing a repeat; the assertion below caught it
# and it is now built by intersecting candidate sets, which cannot.
deg = A.sum(axis=1)
quad = []
for start in np.argsort(-deg)[:40]:
    cur = [int(start)]
    cand = np.nonzero(A[start])[0]
    while cand.size:
        nxt = int(cand[np.argmax(deg[cand])])
        cur.append(nxt)
        cand = cand[A[nxt, cand]]
        cand = cand[~np.isin(cand, cur)]
    if len(cur) > len(quad):
        quad = cur
    if len(quad) >= k:
        break
quad = quad[:k]
for x in range(len(quad)):
    for y in range(x + 1, len(quad)):
        assert A[quad[x], quad[y]], "the clique must be pairwise adjacent"
assert len(set(quad)) == len(quad), "no repeats"
assert len(quad) >= 3, "need at least a triangle to break anything"
print(f"[{name}] clique of {len(quad)} at {quad}, pairwise adjacent in the "
      f"combined graph, so fixing them to {list(range(len(quad)))} is without "
      f"loss of generality -- a factor of "
      f"{np.prod([k - t for t in range(len(quad))])}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in allp:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
for pos, v in enumerate(quad):
    cls.append([1 + v * k + pos])
print(f"[{name}] Sa at five, {tag} classes {sorted(CS)}: {n} pts, {len(allp)} constraints, "
      f"{len(cls)} clauses, colour symmetry broken {np.prod([k - t for t in range(len(quad))])}-fold"
      f"  [{time.time()-t0:.0f}s]", flush=True)
sv = Solver(name=name, bootstrap_with=cls)
ok = sv.solve()
verdict = ("colours -- SATISFIABLE" if ok else
           "*** DOES NOT COLOUR -- THE WEAK PROPERTY HOLDS AT FIVE ***")
print(f"[{name}] -> {verdict}  [{time.time()-t0:.0f}s]", flush=True)
