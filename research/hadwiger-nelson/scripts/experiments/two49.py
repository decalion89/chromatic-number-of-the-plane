"""CDCL on the instance sitting at one conflict.

Sa at five colours with classes 4/9 and 16/9 forbidden -- 2523 constraints on
397 points -- is where the calibrated local search gets closest to the line
anywhere in this work.  Ten runs: best ONE conflict, reached four times, never
zero.  For scale, the known-unsatisfiable calibration plateaus at thirty-seven,
which is a badly unsatisfiable instance; one is a barely unsatisfiable one, or
a satisfiable one the search keeps missing by a single vertex.

Local search cannot tell those apart -- it can only ever prove the
satisfiable side, by landing.  CDCL can prove the other.  So put it on this
exact instance rather than the four-class one, with three solvers, and let it
run.

If it comes back UNSAT then in every 5-colouring of Sa some pair at distance
2/3 or 4/3 is monochromatic: the first weak property at five colours anywhere
in this repository.  If it comes back SAT the local search was missing by one
and that is worth knowing too.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k = 5
CS = {Fr(4, 9), Fr(16, 9)}
name = sys.argv[1]
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
cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
for a, b in sorted(E) + extra:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
print(f"[{name}] Sa at five, {{4/9, 16/9}}: {g.n} pts, {len(E)} unit edges + "
      f"{len(extra)} forbidden pairs, {len(cls)} clauses"
      f"  [{time.time()-t0:.0f}s]", flush=True)
sv = Solver(name=name, bootstrap_with=cls)
ok = sv.solve()
verdict = ("colours -- satisfiable, the local search was missing by one"
           if ok else
           "*** DOES NOT COLOUR -- THE WEAK PROPERTY HOLDS AT FIVE COLOURS ***")
print(f"[{name}] -> {verdict}  [{time.time()-t0:.0f}s]", flush=True)
