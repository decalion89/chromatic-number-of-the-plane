"""Mixed-distance machinery: denesting, conflict rotations, escape counts."""
from fractions import Fraction

from hn.degrey import build_Sa
from hn.field import Field
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.mixed import (_sqrt_in_field, conflict_rotation_set,
                      count_cross_transversals)


def test_rational_square_roots():
    f = Field((3, 5, 7, 11))
    assert _sqrt_in_field(f.rational(4)) == f.rational(2)
    r = _sqrt_in_field(f.rational(3))
    assert r is not None and r * r == f.rational(3)
    assert _sqrt_in_field(f.rational(2)) is None          # sqrt2 not in here
    assert _sqrt_in_field(f.rational(-1)) is None


def test_denesting_two_term_elements():
    """sqrt(a + b sqrt d) = sqrt(x) + sqrt(y) when a^2 - b^2 d is a square.

    The case that mattered: 7/72 + sqrt(33)/72 comes out as
    (sqrt(11) + sqrt(3))/12, and restricting to rational arguments returned
    nothing at all for the core it belongs to.
    """
    f = Field((3, 5, 7, 11))
    v = f.rational(Fraction(7, 72)) + f.sqrt(33) * f.rational(Fraction(1, 72))
    r = _sqrt_in_field(v)
    assert r is not None and r * r == v
    expect = (f.sqrt(11) + f.sqrt(3)) * f.rational(Fraction(1, 12))
    assert r == expect


def test_denesting_refuses_what_is_not_there():
    """A negative a^2 - b^2 d means the element is not totally positive."""
    f = Field((3, 5, 7, 11))
    v = f.rational(Fraction(11, 12)) + f.sqrt(33) * f.rational(Fraction(1, 6))
    assert (Fraction(11, 12) ** 2 - 33 * Fraction(1, 6) ** 2) < 0
    assert _sqrt_in_field(v) is None


def test_conflict_rotations_do_create_conflicts():
    """Each returned rotation puts that pair's image exactly one apart.

    Most pairs return nothing, since the square root the rotation needs is
    usually outside the field -- 52 of 81 at one measured pivot. So this hunts
    for a pair that does resolve rather than assuming the first few will.
    """
    from hn.mixed import conflict_rotations

    g = build_graph(build_Sa())
    pivot = 0
    p = g.vertices[pivot]
    cand = [j for j in range(1, g.n) if j not in g.adj[pivot]][:40]
    found = 0
    for a in cand:
        for b in cand:
            for r in conflict_rotations(g, pivot, a, b):
                assert r.about(p)(g.vertices[a]).is_unit_apart(g.vertices[b])
                found += 1
        if found >= 3:
            break
    assert found, "no nameable conflict rotation anywhere in the sample"


def test_escape_count_is_zero_only_when_blocked():
    from hn.multispindle import cross_blocks

    g = build_graph(build_Sa())
    field = g.vertices[0].x.field
    pivot = 0
    targets = [j for j in range(1, g.n) if j not in g.adj[0]][:3]
    rots = [Rotation(field.rational(1), field.zero())]
    rots += conflict_rotation_set(g, pivot, targets)[:4]
    n = count_cross_transversals(g, pivot, targets, rots, cap=5000)
    assert (n == 0) == cross_blocks(g, pivot, targets, rots)


def test_unit_circle_intersections_are_exactly_one_away():
    """The points at distance 1 from both ends of a short edge."""
    from hn.geometry import eisenstein, origin
    from hn.mixed import unit_circle_intersections

    u, v = origin(), eisenstein(1, 0)          # one apart
    xs = unit_circle_intersections(u, v)
    assert len(xs) == 2
    for x in xs:
        assert x.dist2(u) == 1 and x.dist2(v) == 1
    assert xs[0] != xs[1]


def test_no_intersection_when_the_ends_are_too_far():
    from hn.geometry import eisenstein, origin
    from hn.mixed import unit_circle_intersections

    u, v = origin(), eisenstein(3, 0)          # three apart: no such point
    assert unit_circle_intersections(u, v) == []


def test_deep_holes_really_have_the_degree_claimed():
    """A pivot need not be a vertex; these are the points worth adding."""
    from hn.mixed import deep_holes

    g = build_graph(build_Sa())
    holes = deep_holes(g, min_degree=8, limit=200)
    assert holes
    for deg, x in holes[:5]:
        assert sum(1 for q in g.vertices if x.is_unit_apart(q)) == deg
