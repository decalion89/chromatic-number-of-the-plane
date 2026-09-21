"""Bite the class that is already at the edge.

The by-distance gateway is a far more sensitive instrument than the antipodal
one, and it costs accordingly.  G at five colours with every pair at squared
distance 4/9 forbidden resisted four solvers for over an hour, where the
classes with 6510, 3648, 3216 and 2448 pairs each fell in seconds.  Balls
about the hub place it exactly: 1200 points in 21s, 1400 in 72s, 1500 in
1279s, 1540 in 4969s, 1560 in 6893s -- and all of them COLOUR, with 1548 of
the 1558 pairs forbidden.  So G almost certainly colours too, and the useful
reading is not the verdict but the curve: G sits right at the edge of the
property on this class, and nowhere near it on any other.

D = 4/9 is also doubly usable -- 4D - 1 = 7/9 wants sqrt(7) for the bite and
16D - 1 = 55/9 wants sqrt(55) = sqrt(5) sqrt(11) for the spindle, both already
in de Grey's field -- and G carries a six-point ring of it about 1239 of its
vertices.  Which makes this the best-motivated construction available: the
template's own move, applied to the one class the graph is already straining
on, with no adjunction needed.

Bite it both ways about a centre that has the ring, and ask the class again.
Balls first, because the full instance costs hours and the curve is what
answers: if the union's curve sits left of G's, the bite is tightening the
class, and if it crosses, the property holds.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining, Rotation
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import doubly_usable_ring
from pysat.solvers import Solver

k, D = 5, Fr(4, 9)
t0 = time.time()
print(f"D={D} doubly usable: {doubly_usable_ring(D)}", flush=True)
G = build_G(F, as_graph=False)
# a centre that carries the ring
C = None
for i, p in enumerate(G):
    n = sum(1 for q in G if q.dist2(p) == D)
    if n >= 6:
        C, ci, ring = p, i, n
        break
print(f"centre G[{ci}] = ({C.fx:+.4f}, {C.fy:+.4f}) with {ring} points on "
      f"the D={D} ring  [{time.time()-t0:.0f}s]", flush=True)
b = rotation_joining(D, F)
print(f"bite: cos {b.cos}, sin {b.sin}", flush=True)
f, iv = b.about(C), Rotation(b.cos, -b.sin).about(C)
q = next(p for p in G if p.dist2(C) == D)
print(f"check: a ring point and its image are {q.dist2(f(q))} apart",
      flush=True)

seen, U = set(), []
for p in list(G) + [f(p) for p in G] + [iv(p) for p in G]:
    if p not in seen:
        seen.add(p)
        U.append(p)
print(f"union: {len(U)} points  [{time.time()-t0:.0f}s]", flush=True)

r2 = np.array([float(p.dist2(C)) for p in U])
order = np.argsort(r2)
print("\nthe curve, against G's 21s/72s/1279s/4969s/6893s at "
      "1200/1400/1500/1540/1560:", flush=True)
for m in (600, 900, 1200, 1400, 1600, 1900, 2300, 2800, 3400, len(U)):
    if m > len(U):
        break
    keep = sorted(int(i) for i in order[:m])
    Q = [U[i] for i in keep]
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
    for a, bb in sorted(E) + pr:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + bb * k + c)])
    t1 = time.time()
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    verdict = ("colours" if ok else
               "*** DOES NOT COLOUR -- THE PROPERTY HOLDS ***")
    print(f"  {g.n} pts, {len(E)} edges, {len(pr)} pairs at D={D} -> "
          f"{verdict} in {time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]",
          flush=True)
    if not ok:
        break
print("DONE", flush=True)
