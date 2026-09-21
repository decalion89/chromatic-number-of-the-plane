"""Balls, again, because a witness is smaller than a proof.

Sa at five with all four carrying classes forbidden -- 825 pairs on 397 points
-- has held out against cadical and glucose.  A 397-vertex 5-colouring
instance normally resolves in milliseconds, so it is at the phase transition.

Forbidding is monotone: if a SUBGRAPH does not colour with its share of those
pairs forbidden, neither does Sa, and the subgraph is the witness.  Small balls
answer instantly, so climb them.  Every ball that colours costs nothing; the
first that does not is the answer AND the smallest object carrying a weak
property at five colours anywhere here.

Sa is the dihedral closure about the origin, so balls about the origin respect
its symmetry rather than cutting across it.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k = 5
CARRY = {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}
t0 = time.time()
which = sys.argv[1] if len(sys.argv) > 1 else "Sa"
P0 = build_Sa(F) if which == "Sa" else build_Y(F)
O = Point(F.zero(), F.zero())
r2 = np.array([float(p.dist2(O)) for p in P0])
order = np.argsort(r2)
print(f"{which}: {len(P0)} points about the origin  [{time.time()-t0:.0f}s]",
      flush=True)

sizes = [40, 60, 80, 100, 130, 160, 190, 220, 250, 280, 310, 340, 370,
         len(P0)]
for m in sizes:
    if m > len(P0):
        break
    Q = [P0[int(i)] for i in sorted(order[:m])]
    g = build_graph(Q)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(Q)
    rows = basis.rows(Q)
    dim, D2 = basis.dim, basis.D * basis.D
    pr = []
    for i in range(len(Q) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for t in range(1, dim):
            rat &= sq[:, t] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) in CARRY:
                pr.append((i, j))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E) + pr:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    t1 = time.time()
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    verdict = ("colours" if ok else
               "*** DOES NOT COLOUR -- WITNESS AT FIVE COLOURS ***")
    print(f"  ball of {g.n} pts, {len(E)} edges, {len(pr)} pairs forbidden "
          f"-> {verdict} in {time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]",
          flush=True)
    if not ok:
        break
print("DONE", flush=True)
