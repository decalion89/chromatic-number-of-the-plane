"""Verify the colouring, and note what the local search was telling us.

glucose found a 5-colouring of Sa in which every pair at squared distance 4/9,
16/9, 4 and 16 is bichromatic -- 850 seconds, colour symmetry broken by a
triangle.  So Sa does NOT carry the weak property at five colours on those
four classes.

The two-class statement follows for free and needs no run of its own: the
forbidden set {4/9, 16/9} is a SUBSET of {4/9, 16/9, 4, 16}, so a colouring
avoiding the larger set avoids the smaller one.  Satisfiable there too.

Which retires the evidence that pointed the other way.  Sixty-seven runs of a
calibrated TabuCol -- sixty at 600k moves and seven at three million -- never
reached zero on the two-class instance and touched ONE conflict twenty-four
times.  That looked like a barely unsatisfiable instance.  It was a hard
satisfiable one, and the same thing had already happened once, on G's class
4/9, where a three-hour CDCL standoff also ended in SAT.  A local search that
does not land is not evidence; only landing proves anything, and only in one
direction.

The check here is against the exact geometry, not the solver's word.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, time
sys.path.insert(0, HN_DIR)
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k = 5
CS = {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}
t0 = time.time()
P = build_Sa(F)
g = build_graph(P)
E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
Eset = set(E)
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
        if (i, j) not in Eset and Fr(int(sq[off, 0]), D2) in CS:
            extra.append((i, j))
print(f"Sa: {g.n} pts, {len(E)} unit edges, {len(extra)} forbidden pairs "
      f"over {sorted(CS)}  [{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
for a, b in E + extra:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="glucose4", bootstrap_with=cls)
ok = sv.solve()
print(f"solver says {'SAT' if ok else 'UNSAT'}  [{time.time()-t0:.0f}s]",
      flush=True)
assert ok
mo = sv.get_model()
col = [next(c for c in range(k) if mo[v * k + c] > 0) for v in range(g.n)]
sv.delete()

bad_unit = [(a, b) for a, b in E if col[a] == col[b]]
bad_extra = [(a, b) for a, b in extra if col[a] == col[b]]
wrong_unit = sum(1 for a, b in E if P[a].dist2(P[b]) != 1)
wrong_extra = sum(1 for a, b in extra if P[a].dist2(P[b]) not in CS)
print(f"  vertices coloured: {len(col)} of {g.n}, colours {sorted(set(col))}",
      flush=True)
print(f"  unit edges monochromatic:      {len(bad_unit)}", flush=True)
print(f"  forbidden pairs monochromatic: {len(bad_extra)}", flush=True)
print(f"  edges whose exact distance is not 1:        {wrong_unit}",
      flush=True)
print(f"  pairs whose exact distance is not in the classes: {wrong_extra}",
      flush=True)
ok2 = not (bad_unit or bad_extra or wrong_unit or wrong_extra)
print(f"\n{'VERIFIED' if ok2 else 'DOES NOT CHECK OUT'}: Sa at five colours "
      f"admits a colouring avoiding all four classes.", flush=True)
print("So Sa carries NO weak property at five colours, on those classes "
      "singly or together; and by monotonicity the {4/9, 16/9} statement "
      "is satisfiable too.", flush=True)
print("DONE", flush=True)
