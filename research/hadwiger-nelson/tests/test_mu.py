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


# --- positive controls -------------------------------------------------
#
# Everything above shows mu coming back at its floor.  That is worthless as
# evidence unless the instrument can also come back at the CEILING when the
# structure is really there, so here it is made to.
#
# The Moser spindle is 4-chromatic and 4-vertex-critical, so deleting any
# vertex p leaves a 3-colourable graph in which N(p) must carry all three
# colours -- otherwise p could be put back and the spindle would be
# 3-chromatic.  That is mu_3 = 3 = k: a blocked point, at three colours,
# with a known answer.


def _moser():
    """Two rhombi hinged at the origin, their far apexes a unit apart.

    A rhombus is two unit triangles sharing an edge: (0,0), (sqrt3/2, +-1/2)
    and (sqrt3, 0), with the apexes sqrt3 apart.  Two of them, rotated against
    each other by the angle with cos 5/6 and sin sqrt11/6, put the far apexes
    at 6(1 - 5/6) = 1 -- which is the whole point of the construction and the
    reason sqrt11 is in de Grey's field at all.
    """
    f = Field((3, 11))
    half = f.rational(Fr(1, 2))
    rt3 = f.sqrt(3)
    base = [Point(f.zero(), f.zero()),
            Point(rt3 * half, half),
            Point(rt3 * half, -half),
            Point(rt3, f.zero())]
    spin = Rotation(f.rational(Fr(5, 6)), f.sqrt(11) * f.rational(Fr(1, 6)))
    pts, seen = [], set()
    for p in base + [spin(q) for q in base]:
        if p not in seen:
            seen.add(p); pts.append(p)
    return f, pts


def test_the_moser_spindle_is_what_it_should_be():
    f, pts = _moser()
    g = build_graph(pts)
    assert g.n == 7
    assert sum(len(a) for a in g.adj) // 2 == 11
    one = f.rational(Fr(1))
    assert (pts[3] - pts[6]).norm2() == one, "the two apexes must be adjacent"


@pytest.mark.parametrize("drop", range(7))
def test_every_deleted_vertex_of_the_spindle_is_blocked_at_three(drop):
    """mu reaches the ceiling, on all seven vertices.

    The spindle is 4-vertex-critical, so this must hold for every choice of p,
    and if the instrument ever reported less than 3 it would be under-counting
    -- which is the failure mode that matters, since every result at five
    colours here is a claim that mu is SMALL.
    """
    _, pts = _moser()
    p = pts[drop]
    rest = build_graph([q for i, q in enumerate(pts) if i != drop])
    nb = neighbourhood(rest, p)
    assert len(nb) >= 3
    with MuSolver(rest, k=3) as ms:
        assert ms.colourable, "the spindle minus a vertex is 3-colourable"
        assert ms.mu(nb) == 3
        assert ms.at_most_two(nb) is False


def test_the_spindle_is_not_blocked_at_four():
    """And mu must not over-report: at four colours the same point is placeable
    again, since the spindle itself is 4-chromatic."""
    _, pts = _moser()
    p = pts[3]
    rest = build_graph([q for i, q in enumerate(pts) if i != 3])
    nb = neighbourhood(rest, p)
    with MuSolver(rest, k=4) as ms:
        assert ms.mu(nb) < 4
