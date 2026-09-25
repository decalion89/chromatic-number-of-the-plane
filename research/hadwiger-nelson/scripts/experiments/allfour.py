"""Skip to the most constrained case: all four carrying classes at once.

Sa at five with {4/9, 16/9} forbidden has been running fourteen minutes on a
397-point graph, where 4/9 alone took ten seconds.  A factor of eighty on an
instance that small means it is at the phase transition, not merely slow.

More constraints usually make an UNSAT proof EASIER, not harder -- the search
space collapses faster.  So rather than climb through the pairs and triples,
ask the most constrained question directly: forbid all four classes that Sa
carries at four colours -- 4/9, 16/9, 4 and 16, which is 393 + 156 + 273 + 3
pairs -- and ask for a 5-colouring.

If it comes back UNSAT then in every 5-colouring of Sa some pair at distance
2/3, 4/3, 2 or 4 is monochromatic.  That is four distances rather than one, so
the spindle that follows is a multispindle over four rotation angles rather
than one -- awkward, not impossible -- and it would be the first weak property
at five colours anywhere in this repository.

Two solvers, since the answer matters and the instance is small.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

which, name = sys.argv[1], sys.argv[2]
k = 5
CARRY = {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}
t0 = time.time()
P = build_Sa(F) if which == "Sa" else build_Y(F)
g = build_graph(P)
E = set((min(a, b), max(a, b)) for a, b in g.edges())
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
byd = defaultdict(list)
for i in range(len(P) - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in E:
            byd[Fr(int(sq[off, 0]), D2)].append((i, j))
pr = [p for d in CARRY for p in byd.get(d, [])]
sizes = {str(d): len(byd.get(d, [])) for d in sorted(CARRY)}
cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
for a, b in sorted(E) + pr:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
print(f"[{name}] {which}: {g.n} pts, {len(E)} edges, classes {sizes}, "
      f"{len(pr)} pairs forbidden  [{time.time()-t0:.0f}s]", flush=True)
sv = Solver(name=name, bootstrap_with=cls)
ok = sv.solve()
verdict = ("colours (satisfiable)" if ok else
           "*** DOES NOT COLOUR -- THE WEAK PROPERTY HOLDS AT FIVE ***")
print(f"[{name}] {which} -> {verdict}  [{time.time()-t0:.0f}s]", flush=True)
