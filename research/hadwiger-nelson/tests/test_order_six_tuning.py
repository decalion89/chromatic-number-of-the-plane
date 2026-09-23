"""Tuning the ANGLE instead of the distance closes the chain into a group.

The composition tau(p) = q + R_phi (p - v) is a rotation, since its linear part
is one.  It is therefore a rotation by phi about the point w it fixes, and
because tau(v) = q the points v and q sit on a circle about w separated by phi:

        |v - q|^2 = 2 R^2 (1 - cos phi),      R = |w - v|

Every earlier use chose phi to place the composite DISTANCE.  Choosing it to
fix the ORDER gives something the carrier did not have.  At phi = 60 degrees --
which needs only sqrt3, already in every field here -- tau^6 = id, and

  * the union of the six copies tau^j(H) is C6-invariant about w, a centre that
    is not a vertex of the carrier and a group the carrier did not have;
  * the forcing chains all the way round.  tau^j(H) forces the pair
    (tau^j v, tau^j q), and tau^j v = tau^{j-1} q, so consecutive points of the
    ORBIT of v are forced equal and, by transitivity, all six take one colour;
  * R^2 = D^2 / (2(1 - cos 60)) = D^2, so with D^2 = 64/9 the orbit is six
    points on a circle of radius 8/3 and its fifteen squared distances are
    64/9, 64/3 and 256/9 -- six, six and three of them.

That last line is why it matters.  A forced class at three different radii is
exactly what the free-angle spindle needs, and it had turned up once, by luck,
on one 3025-point carrier.  This produces it from ANY forced pair at 64/9.

Of the three radii only 64/3 with 256/9 admits the free angle:
|sqrt(64/3) - sqrt(256/9)| = 0.715 <= 1, so the circles meet; 64/9 against
either gives 1.952 and they do not.
"""
from __future__ import annotations

import json
import os
from collections import Counter
from fractions import Fraction as Fr

import pytest

from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

HERE = os.path.dirname(__file__)
FIELD = Field((3, 11, 247))
D2 = Fr(64, 9)
COS = FIELD.rational(Fr(1, 2))
SIN = FIELD.sqrt(3) * FIELD.rational(Fr(1, 2))


@pytest.fixture(scope="module")
def stored():
    with open(os.path.join(HERE, os.pardir, "data",
                           "order6_carrier.json")) as fh:
        d = json.load(fh)
    pts = [Point(FIELD.element([Fr(a, b) for a, b in x]),
                 FIELD.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


def _tau(v, q):
    r = Rotation(COS, SIN)

    def f(p):
        d = Point(p.x - v.x, p.y - v.y)
        z = r(d)
        return Point(q.x + z.x, q.y + z.y)
    return f


def test_the_angle_is_sixty_and_the_order_is_six(stored):
    d, pts = stored
    orbit = [pts[i] for i in d["forced_orbit"]]
    assert len(set(orbit)) == 6
    v, q = orbit[0], orbit[1]
    assert (v - q).norm2() == FIELD.rational(D2)
    tau = _tau(v, q)
    z = v
    for _ in range(6):
        z = tau(z)
    assert z == v, "tau^6 must fix v exactly"
    # and no smaller power does
    z = v
    for k in range(1, 6):
        z = tau(z)
        assert z != v, f"tau^{k} must not fix v"


def test_the_fixed_point_sits_at_the_forced_distance(stored):
    """R^2 = D^2/(2(1 - cos phi)) = D^2 at phi = 60, so the orbit radius IS the
    forced distance.

    w is recomputed from scratch.  tau(w) = w means (I - R)w = q - Rv, and

        I - R = [[1-c, s], [-s, 1-c]],   det = (1-c)^2 + s^2 = 2(1-c)

    which is exactly 1 at c = 1/2, so the inverse is [[1-c, -s], [s, 1-c]] --
    that is [[1/2, -sqrt3/2], [sqrt3/2, 1/2]], the rotation by +60 degrees, NOT
    by -60.  Using the inverse rotation instead put w off the circle.
    """
    d, pts = stored
    orbit = [pts[i] for i in d["forced_orbit"]]
    v, q = orbit[0], orbit[1]
    fwd = Rotation(COS, SIN)
    rv = fwd(v)
    w = fwd(Point(q.x - rv.x, q.y - rv.y))
    assert (w - v).norm2() == FIELD.rational(D2)
    for p in orbit:
        assert (w - p).norm2() == FIELD.rational(D2), "orbit is not a circle"


def test_the_orbit_realises_exactly_three_squared_distances(stored):
    d, pts = stored
    orbit = [pts[i] for i in d["forced_orbit"]]
    c = Counter(str((orbit[i] - orbit[j]).norm2())
                for i in range(6) for j in range(i + 1, 6))
    assert dict(c) == {"64/9": 6, "64/3": 6, "256/9": 3}
    assert sorted(c) == sorted(d["orbit_d2"])


def test_only_one_of_the_three_pairs_admits_a_free_angle():
    """|sqrt(r1) - sqrt(r2)| <= 1 is the condition for the two circles to meet,
    which is what the free-angle rotation needs."""
    import math
    rs = [Fr(64, 9), Fr(64, 3), Fr(256, 9)]
    ok = [(a, b) for i, a in enumerate(rs) for b in rs[i + 1:]
          if abs(math.sqrt(a) - math.sqrt(b)) <= 1]
    assert ok == [(Fr(64, 3), Fr(256, 9))]


def test_the_union_is_invariant_and_matches_its_counts(stored):
    d, pts = stored
    g = build_graph(pts)
    assert g.n == d["n"]
    assert sum(len(a) for a in g.adj) // 2 == d["m"]
    orbit = [pts[i] for i in d["forced_orbit"]]
    tau = _tau(orbit[0], orbit[1])
    seen = set(g.vertices)
    assert all(tau(p) in seen for p in g.vertices), "not tau-invariant"


@pytest.mark.parametrize("i,j", [(0, 2), (0, 3)])
def test_the_orbit_is_forced_equal(stored, i, j):
    """Two of the fifteen pairs, one at 64/3 and one at 256/9, asked directly.

    The CNF carries no symmetry breaking, so "x takes colour 0 and y colour 1
    is impossible" does mean they agree in every 4-colouring.
    """
    d, pts = stored
    g = build_graph(pts)
    K = 4
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(g.n)]
    for v in range(g.n):
        for a in range(K):
            for b in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    a, b = d["forced_orbit"][i], d["forced_orbit"][j]
    s = Solver(name="cd19", bootstrap_with=cnf)
    assert s.solve(), "the carrier must be 4-colourable"
    forced = not s.solve(assumptions=[X(a, 0), X(b, 1)])
    s.delete()
    assert forced
