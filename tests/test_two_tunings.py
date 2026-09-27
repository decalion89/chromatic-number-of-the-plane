"""Both tunings at once: the group and the chromatic number, together.

Tuning the DISTANCE moves a forced pair anywhere in [0, 2D], at the cost of one
radical.  Tuning the ANGLE to 60 degrees makes the composition have order six,
so the union of its copies is C6-invariant and the whole orbit of the pivot
takes one colour, on a circle of radius R with R^2 = D^2.

Applied in that order they collide usefully.  Move the pair to squared distance
1/3 first -- cos phi1 = -125/128, sin phi1 = sqrt(759)/128 with 759 = 3*11*23 --
and then tune the angle on the NEW pair.  The orbit now lies on a circle of
radius 1/sqrt3 with its six points 60 degrees apart, so the pairs two apart are
at distance

        2 R sin 60 = 2 * (1/sqrt3) * (sqrt3/2) = 1

exactly.  Six of the fifteen orbit pairs are EDGES, and every one of the
fifteen is forced monochromatic.  The union of the six copies therefore has no
4-colouring, with no spindle anywhere in the construction, and it carries C6
about a centre belonging to none of the copies.

The distance multiset is the whole proof in one line: {1/3: 6, 1: 6, 4/3: 3}.
"""
from __future__ import annotations

import json
import os
from collections import Counter
from fractions import Fraction as Fr

import pytest

from hn.coloring import is_k_colorable
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

HERE = os.path.dirname(__file__)
FIELD = Field((3, 11, 23))


@pytest.fixture(scope="module")
def stored():
    with open(os.path.join(HERE, os.pardir, "data", "five_twotune.json")) as fh:
        d = json.load(fh)
    pts = [Point(FIELD.element([Fr(a, b) for a, b in x]),
                 FIELD.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


def test_both_tuning_angles_are_exact_rotations():
    c1 = FIELD.rational(Fr(-125, 128))
    s1 = FIELD.sqrt(759) * FIELD.rational(Fr(1, 128))
    assert c1 * c1 + s1 * s1 == FIELD.rational(Fr(1))
    # and it is the angle that lands 64/9 on 1/3
    assert c1 == FIELD.rational(Fr(1, 3) / (2 * Fr(64, 9)) - 1)
    c2 = FIELD.rational(Fr(1, 2))
    s2 = FIELD.sqrt(3) * FIELD.rational(Fr(1, 2))
    assert c2 * c2 + s2 * s2 == FIELD.rational(Fr(1))


def test_the_distance_tuning_lands_on_one_third():
    """|u + R u|^2 = 2 D^2 (1 + cos phi), checked on a vector of length D."""
    c1 = FIELD.rational(Fr(-125, 128))
    s1 = FIELD.sqrt(759) * FIELD.rational(Fr(1, 128))
    u = Point(FIELD.rational(Fr(8, 3)), FIELD.zero())
    assert u.norm2() == FIELD.rational(Fr(64, 9))
    r = Rotation(c1, s1)(u)
    assert Point(u.x + r.x, u.y + r.y).norm2() == FIELD.rational(Fr(1, 3))


def test_the_orbit_puts_six_pairs_at_distance_exactly_one(stored):
    d, pts = stored
    orbit = [pts[i] for i in d["orbit"]]
    assert len(set(orbit)) == 6
    c = Counter(str((orbit[i] - orbit[j]).norm2())
                for i in range(6) for j in range(i + 1, 6))
    assert dict(c) == {"1/3": 6, "1": 6, "4/3": 3}
    assert dict(c) == d["orbit_d2"]
    # the unit pairs are the ones two apart on the circle
    for i in range(6):
        assert (orbit[i] - orbit[(i + 2) % 6]).norm2() == FIELD.rational(Fr(1))


def test_the_orbit_lies_on_a_circle_of_radius_one_over_root_three(stored):
    """R^2 = D^2 at phi = 60, and D^2 is now 1/3 -- which is exactly the radius
    whose inscribed equilateral triangle has side 1."""
    d, pts = stored
    orbit = [pts[i] for i in d["orbit"]]
    c2 = FIELD.rational(Fr(1, 2))
    s2 = FIELD.sqrt(3) * FIELD.rational(Fr(1, 2))
    fwd = Rotation(c2, s2)
    v, q = orbit[0], orbit[1]
    rv = fwd(v)
    w = fwd(Point(q.x - rv.x, q.y - rv.y))          # (I - R)^-1 = R at 60
    for p in orbit:
        assert (w - p).norm2() == FIELD.rational(Fr(1, 3))


def test_the_union_is_c6_invariant_and_matches_its_counts(stored):
    d, pts = stored
    g = build_graph(pts)
    assert g.n == d["n"]
    assert sum(len(a) for a in g.adj) // 2 == d["m"]
    orbit = [pts[i] for i in d["orbit"]]
    v, q = orbit[0], orbit[1]
    r = Rotation(FIELD.rational(Fr(1, 2)), FIELD.sqrt(3) * FIELD.rational(Fr(1, 2)))

    def tau(p):
        z = r(Point(p.x - v.x, p.y - v.y))
        return Point(q.x + z.x, q.y + z.y)
    seen = set(g.vertices)
    assert all(tau(p) in seen for p in g.vertices)


def test_it_refuses_four(stored):
    """The claim itself.  Slow -- the refutation is over 12 240 vertices and the
    solver has to rediscover the forcing that makes the six unit pairs
    monochromatic."""
    _, pts = stored
    g = build_graph(pts)
    answer, colouring = is_k_colorable(g, 4, timeout=4 * 3600)
    assert answer is False
    assert colouring is None
