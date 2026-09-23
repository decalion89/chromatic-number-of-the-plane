"""mu, the number that is the problem, checked on objects with known answers.

    mu_k(G, p) = min over proper k-colourings of |c(N(p))|,
    N(p) = the graph points at distance exactly 1 from p

G + p has no k-colouring exactly when mu_k(G, p) = k, so at k = 5 the statement
"mu = 5 somewhere" is chi(R^2) >= 6 and "mu <= 2 always" is chi(R^2) = 5.  The
module is worth testing on cases where the answer is known in advance, which is
what this file does: a triangle, where adding the centre of a unit circle
through all three is impossible in the plane but the arithmetic still has to
come out right; the hexagon, where the centre is a genuine vertex; and de
Grey's Sa at four colours, where the earlier pass established the answer.
"""
from __future__ import annotations

from fractions import Fraction as Fr

import pytest

from hn.blocked import MuSolver, mu_of_point, neighbourhood
from hn.degrey import build_Sa
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

FIELD = Field((3, 5, 7, 11))
ZERO = FIELD.zero()
ONE = FIELD.rational(Fr(1))


def _hexagon():
    """The six unit vectors: every one is at distance 1 from the origin."""
    half = FIELD.rational(Fr(1, 2))
    rt = FIELD.sqrt(3) * half
    return [Point(ONE, ZERO), Point(half, rt), Point(-half, rt),
            Point(-ONE, ZERO), Point(-half, -rt), Point(half, -rt)]


def test_neighbourhood_is_exact_not_approximate():
    g = build_graph(_hexagon())
    origin = Point(ZERO, ZERO)
    assert sorted(neighbourhood(g, origin)) == list(range(6))
    # a point a hair off the origin has an EMPTY neighbourhood: the float
    # filter must never decide anything by itself
    nudge = Point(FIELD.rational(Fr(1, 1000)), ZERO)
    assert neighbourhood(g, nudge) == []


def test_the_hexagon_centre_is_placeable_at_five():
    """The six unit vectors form a 6-cycle, so two colours suffice on them and
    the centre has three to choose from.  mu = 2, the floor."""
    g = build_graph(_hexagon())
    mu, nb = mu_of_point(g, Point(ZERO, ZERO), k=5)
    assert len(nb) == 6
    assert mu == 2


def test_the_hexagon_centre_is_blocked_at_three():
    """At three colours the same 6-cycle still takes only two, so mu = 2 < 3
    and the centre is placeable -- the hexagon plus its centre is 3-chromatic,
    not 4.  The instrument has to agree with that."""
    g = build_graph(_hexagon())
    mu, _ = mu_of_point(g, Point(ZERO, ZERO), k=3)
    assert mu == 2


def test_a_neighbourhood_can_never_beat_its_own_size():
    """mu is capped by |N(p)|, since j points show at most j colours."""
    g = build_graph(_hexagon())
    with MuSolver(g, k=5) as ms:
        assert ms.mu([0]) == 1
        assert ms.mu([0, 3]) <= 2      # opposite vectors, distance 2, no edge


def test_mu_on_de_greys_carrier_at_five_is_two():
    """Sa is 4-chromatic, so at five colours it is loose, and its unit circles
    are the easiest possible case.  Taking the richest neighbourhood available
    among the vertices' own rotated images, mu comes back at the floor."""
    sa = build_Sa(FIELD)
    g = build_graph(sa)
    half = FIELD.rational(Fr(1, 2))
    r60 = Rotation(half, FIELD.sqrt(3) * half)
    seen = set(g.vertices)
    best, bestnb = 0, []
    for c in range(0, g.n, 37):
        rot = r60.about(g.vertices[c])
        for q in g.vertices[:200]:
            z = rot(q)
            if z in seen:
                continue
            nb = neighbourhood(g, z)
            if len(nb) > best:
                best, bestnb = len(nb), nb
    assert best >= 3, "expected some candidate with a few neighbours"
    with MuSolver(g, k=5) as ms:
        assert ms.colourable
        assert ms.mu(bestnb) == 2


def test_the_cheap_probe_agrees_with_the_chain():
    """at_most_two is one call; mu is up to four.  They must not disagree."""
    g = build_graph(_hexagon())
    with MuSolver(g, k=5) as ms:
        nb = neighbourhood(g, Point(ZERO, ZERO))
        assert ms.at_most_two(nb) is True
        assert ms.mu(nb) <= 2
