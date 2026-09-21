"""The pair de Grey deletes sits at the family's distinguished distance.

Three measurements put D = 4/9 -- distance 2/3 -- apart from every other
closable class in this family.  Sa CARRIES the weak property there at four
colours; at five, Sa's 4/9 call is the only one of fourteen that takes ten
seconds instead of none; and G's 4/9 call at five is the only one of its
classes to resist four solvers for an hour.

build_Y deletes exactly two points from Sa u Sb: (1/3, 0) and (-1/3, 0).
Their squared distance is 4/9.  They are the antipodal pair of the D = 1/9
ring about the origin, and one of the 393 pairs in the class.

That is a coincidence worth pinning rather than interpreting -- his reason for
the deletion is his.  What can be measured is what it costs: whether Y, after
the deletion, still carries the property on that class at four colours.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Sb, build_Y
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

t0 = time.time()


def pairs_at(P, D):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    out = []
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) == D:
                out.append((i, j))
    return g.n, E, out


def ask(name, P, D, k):
    n, E, pr = pairs_at(P, D)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E) + pr:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    t1 = time.time()
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    print(f"  {name} at {k}, D={D}: {n} pts, {len(pr)} pairs -> "
          f"{'colours' if ok else '*** DOES NOT COLOUR -- CARRIES IT ***'}"
          f" in {time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]", flush=True)
    return not ok


D = Fr(4, 9)
SaSb = []
seen = set()
for p in build_Sa(F) + build_Sb(F):
    if p not in seen:
        seen.add(p)
        SaSb.append(p)
print(f"Sa u Sb has {len(SaSb)} points, Y has {len(build_Y(F))}: the two "
      f"deleted are at squared distance "
      f"{Point(F.rational(Fr(1,3)), F.zero()).dist2(Point(F.rational(Fr(-1,3)), F.zero()))}",
      flush=True)
ask("Sa", build_Sa(F), D, 4)
ask("Sa u Sb (before the deletion)", SaSb, D, 4)
ask("Y  (after the deletion)", build_Y(F), D, 4)
ask("Y", build_Y(F), D, 5)
print("DONE", flush=True)
