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


def test_reflections_fix_the_pivot_and_preserve_distance():
    """So the spindle lemma applies to them exactly as to rotations."""
    from hn.geometry import origin
    from hn.mixed import Reflection

    g = build_graph(build_Sa())
    field = g.vertices[0].x.field
    p, q = g.vertices[0], g.vertices[5]
    r = Reflection(field.rational(1), field.zero())
    about = r.about(p)
    assert about(p) == p                      # fixes the pivot
    assert about(q).dist2(p) == q.dist2(p)    # stays on its circle
    assert about(about(q)) == q               # an involution
    assert about(q) != q                      # and not the identity


def test_conflict_reflections_share_the_rotations_discriminant():
    """The square root is the same, so reflections cost nothing to name.

    <sigma(a), b> expands with P = ax bx - ay by and Q = ay bx + ax by, and
    P^2 + Q^2 is still |a|^2 |b|^2 while R is unchanged -- so wherever a
    conflict rotation can be named, a conflict reflection can be too.
    """
    from hn.mixed import conflict_reflections, conflict_rotations

    g = build_graph(build_Sa())
    p = g.vertices[0]
    cand = [j for j in range(1, g.n) if j not in g.adj[0]][:30]
    both = 0
    for a in cand:
        for b in cand:
            rots = conflict_rotations(g, 0, a, b)
            refs = conflict_reflections(g, 0, a, b)
            if rots:
                assert refs, "a rotation resolved but its reflection did not"
                both += 1
            for r in refs:
                assert r.about(p)(g.vertices[a]).is_unit_apart(g.vertices[b])
        if both >= 3:
            break
    assert both


def test_reflections_cut_the_escapes():
    """Measured on the core of three: 432 escapes become 135."""
    from hn.geometry import Rotation
    from hn.mixed import (conflict_isometries, conflict_rotation_set,
                          count_cross_transversals)

    g = build_graph(build_Sa())
    field = g.vertices[0].x.field
    ident = Rotation(field.rational(1), field.zero())
    pivot = 0
    targets = [j for j in range(1, g.n) if j not in g.adj[0]][:3]
    rots = [ident] + conflict_rotation_set(g, pivot, targets)
    both = [ident] + conflict_isometries(g, pivot, targets)
    assert len(both) >= len(rots)
    a = count_cross_transversals(g, pivot, targets, rots, cap=50000)
    b = count_cross_transversals(g, pivot, targets, both, cap=50000)
    assert b <= a


def test_forbidden_patterns_generalise_forcing():
    """A forced disjunction is one forbidden partition; there are others.

    "Some target carries the pivot's colour" says exactly that the partition
    putting the pivot in a block of its own is unrealisable. Nothing restricts
    the question to that shape, and on the Moser spindle at k=4 a three-point
    set already has two of its five partitions forbidden.
    """
    from hn.geometry import SPINDLE, eisenstein, origin
    from hn.mixed import _bell, forbidden_patterns, pattern_pressure

    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    g = build_graph(rh + [SPINDLE(p) for p in rh])
    assert _bell(3) == 5 and _bell(4) == 15
    forb = forbidden_patterns(g, 4, [0, 3, 5])
    assert 0 < len(forb) < _bell(3)
    for part in forb:
        assert sorted(i for b in part for i in b) == [0, 1, 2]
    assert 0 < pattern_pressure(g, 4, [0, 3, 5]) < 1


def test_forbidden_patterns_are_vacuous_without_a_colouring():
    """The spindle is 4-chromatic, so at k=3 every partition is 'forbidden'."""
    from hn.geometry import SPINDLE, eisenstein, origin
    from hn.mixed import _bell, forbidden_patterns

    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    g = build_graph(rh + [SPINDLE(p) for p in rh])
    assert len(forbidden_patterns(g, 3, [0, 3, 5])) == _bell(3)
