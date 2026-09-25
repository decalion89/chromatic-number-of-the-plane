"""D = 4/9 on balls around the hub: a faster answer, and a smaller witness.

The full call -- G at five with all 1558 pairs at squared distance 4/9
forbidden -- has held out against four solvers for the best part of an hour,
where the classes with 6510, 3648, 3216 and 2448 pairs each fell in seconds.
That asymmetry is the interesting part: the smallest of the five instances is
the hard one, and D = 4/9 is one of the three classes Sa carries at four.

Forbidding is monotone in the graph: if a SUBGRAPH does not colour with its
D = 4/9 pairs forbidden, neither does G, and the subgraph is a smaller witness
into the bargain.  So take balls of growing radius about the hub G[0], where
G's ring structure is concentrated, and ask each.  A small ball answers in
seconds, and if one of them says no, it says no for G as well and names a
witness a few hundred points across rather than 1581.

If they all colour, that is not the answer for G -- monotonicity runs one way
-- but it bounds where any witness could live.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k, D = 5, Fr(4, 9)
t0 = time.time()
P = build_G(F, as_graph=False)
C = P[0]
r2 = np.array([float(p.dist2(C)) for p in P])
order = np.argsort(r2)
print(f"G: {len(P)} pts about the hub; radii^2 from {r2[order[0]]:.2f} to "
      f"{r2[order[-1]]:.2f}  [{time.time()-t0:.0f}s]", flush=True)

for m in (1250, 1300, 1350, 1400, 1450, 1500, 1540, 1560, 1570, 1575, 1578, 1580, 1581):
    keep = sorted(int(i) for i in order[:m])
    Q = [P[i] for i in keep]
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
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) == D:
                pr.append((i, j))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E) + pr:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    print(f"  ball of {g.n} pts ({len(E)} edges, {len(pr)} pairs at D={D} "
          f"forbidden) -> "
          f"{'colours' if ok else '*** DOES NOT COLOUR -- WITNESS ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        break
print("DONE", flush=True)
