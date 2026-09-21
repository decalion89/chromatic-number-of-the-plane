"""Race the one hard call: G at five, every pair at squared distance 4/9
forbidden.

The by-distance scan cleared 6510 pairs at D=1/3, 3648 at 7/3, 3216 at 3 and
2448 at 5/9 in seconds each, then stopped dead on D=4/9 with only 1558 pairs.
Twenty-five minutes on the smallest instance of the five is a signal in
itself: either it is UNSAT, which is the property, or it is a hard satisfiable
instance.  D = 4/9 is also one of the three classes Sa carries at four
colours, which is why it is worth not waiting on one solver's luck.

Three solvers, three searches.  A satisfiable instance usually falls quickly
to one of them; agreement on UNSAT is the answer.
"""
import sys, time, os
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k, D = 5, Fr(4, 9)
name = sys.argv[1]
t0 = time.time()
P = build_G(F, as_graph=False)
g = build_graph(P)
E = set((min(a, b), max(a, b)) for a, b in g.edges())
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
pr = []
for i in range(len(P) - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in E and Fr(int(sq[off, 0]), D2) == D:
            pr.append((i, j))
cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
for a, b in sorted(E) + pr:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
print(f"[{name}] {g.n} pts, {len(E)} edges, {len(pr)} pairs at D={D} "
      f"forbidden, {len(cls)} clauses  [{time.time()-t0:.0f}s]", flush=True)
sv = Solver(name=name, bootstrap_with=cls)
ok = sv.solve()
verdict = ("COLOURS (satisfiable)" if ok else
           "*** DOES NOT COLOUR: SOME PAIR AT DISTANCE 2/3 IS ALWAYS "
           "MONOCHROMATIC ***")
print(f"[{name}] -> {verdict}  [{time.time()-t0:.0f}s]", flush=True)
