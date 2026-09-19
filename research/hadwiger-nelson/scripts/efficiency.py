"""Degree efficiency: what a graph actually collects of the density on offer.

Blocking and folding are both monotone in the step set, so they do not trade
off directly -- but density per VERTEX does, and this is the statistic that
shows it.  A module of rank 2 is discrete, so a patch of it is closed: nearly
every point has all its neighbours present.  A module of rank 3 or more is
DENSE in the plane, so no bounded region contains a closed piece, every ball
is boundary, and most of the available steps land outside it.

Write efficiency = average degree / (unit vectors available in the module).
It is the fraction of the possible density a finite construction realises, and
it falls off a cliff exactly where blocking becomes possible.
"""
import sys
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.graph import build_graph
from hn.homcol import edge_vectors, has_homomorphism
from hn.field import Field
from hn.geometry import Point, Rotation

FLD = Field((3, 11))


def tri_lattice(r=5):
    h = FLD.rational(Fraction(1, 2))
    u = Point(h, FLD.sqrt(3) * h)
    one = Point(FLD.one(), FLD.zero())
    return [one.scaled(a) + u.scaled(b)
            for a in range(-r, r + 1) for b in range(-r, r + 1)]


rows = []
for name, g in (("triangular patch", build_graph(tri_lattice())),
                ("de Grey S", build_graph(degrey.build_S())),
                ("de Grey Sa", build_graph(degrey.build_Sa())),
                ("de Grey Y", build_graph(degrey.build_Y())),
                ("de Grey G", degrey.build_G())):
    ev = edge_vectors(g)
    # edge_vectors gives one orientation per edge, and whether both signs turn
    # up depends on how the vertices happen to be indexed -- the triangular
    # lattice came out with 3 where it has 6 neighbours, giving an efficiency
    # over 100%.  Symmetrise, so "steps" is always the number of unit vectors
    # a point could use.
    steps = len({v for d in ev for v in (d, tuple(-x for x in d))})
    avg = 2 * g.m / g.n
    blocked = has_homomorphism(ev, 5)[0] is None
    rows.append((name, g.n, g.m, steps, avg, avg / steps, blocked))

# the blocked balls, from their recorded measurements
rows.append(("Q(zeta_7) ball", 15313, 30276, 174, 2 * 30276 / 15313,
             (2 * 30276 / 15313) / 174, True))
rows.append(("Q(zeta_21) ball", 16015, 32722, 178, 2 * 32722 / 16015,
             (2 * 32722 / 16015) / 178, True))

print(f"{'graph':18} {'n':>6} {'m':>6} {'steps':>6} {'avg deg':>8} "
      f"{'efficiency':>11}  blocked")
print("-" * 66)
for name, n, m, s, avg, eff, b in rows:
    print(f"{name:18} {n:6} {m:6} {s:6} {avg:8.2f} {eff:10.1%}  "
          f"{'YES' if b else 'no'}")
