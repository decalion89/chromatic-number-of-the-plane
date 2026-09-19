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


def test_two_target_blocking_is_two_sat():
    """With two targets each copy makes a binary choice, so escapes are 2-SAT.

    A conflict between two chosen images is a forbidden pair, an escape is a
    satisfying assignment, and the copies block exactly when the instance is
    unsatisfiable. That decides in linear time what cross_blocks searches in
    2^m, so a core of two can be thrown against hundreds of copies at once.
    """
    import random

    from hn.geometry import Rotation
    from hn.mixed import blocks_two_targets, conflict_isometries
    from hn.multispindle import cross_blocks

    g = build_graph(build_Sa())
    field = g.vertices[0].x.field
    ident = Rotation(field.rational(1), field.zero())
    rnd = random.Random(3)
    pool = [j for j in range(1, g.n) if j not in g.adj[0]]
    checked = 0
    for _ in range(40):
        a, b = rnd.sample(pool, 2)
        fam = [ident] + conflict_isometries(g, 0, [a, b])[:6]
        if len(fam) < 2:
            continue
        assert blocks_two_targets(g, 0, [a, b], fam) == cross_blocks(
            g, 0, [a, b], fam)
        checked += 1
    assert checked >= 5


def test_two_target_test_refuses_other_sizes():
    from hn.geometry import Rotation
    from hn.mixed import blocks_two_targets

    g = build_graph(build_Sa())
    field = g.vertices[0].x.field
    ident = Rotation(field.rational(1), field.zero())
    try:
        blocks_two_targets(g, 0, [1, 2, 3], [ident])
    except ValueError:
        return
    raise AssertionError("should refuse a core that is not of size two")


def test_two_orbit_construction_blocks_any_pair_with_a_third_leg():
    """Six copies close any forced pair with one target on d^2 = 1/3.

    A theorem, not a search. On the circle of radius 1/sqrt(3) two points are
    adjacent exactly when 120 degrees apart, so each orbit's three images of a
    form a unit triangle and at most one copy per orbit may choose a --
    leaving at least two choosing b. Two subsets of size two or more of a
    three-element set intersect, so some j has both rho^j and sigma rho^j
    choosing b, and those images differ by sigma, the unit-chord angle of
    their own circle. They are one apart, both carry the pivot's colour, and
    that is the contradiction.

    Checked against both blocking tests, the linear one and the exponential
    one, on every pair the graph offers.
    """
    from fractions import Fraction

    from hn.mixed import blocks_two_targets, two_orbit_block
    from hn.multispindle import cross_blocks

    g = build_graph(build_Sa())
    third = Fraction(1, 3)
    tested = 0
    for bp in range(g.n):
        pv = g.vertices[bp]
        legs = [j for j in range(g.n) if j != bp and j not in g.adj[bp]
                and pv.dist2(g.vertices[j]).is_rational()
                and pv.dist2(g.vertices[j]).c[0] == third]
        if not legs:
            continue
        others = [j for j in range(g.n) if j != bp and j not in g.adj[bp]
                  and j not in legs
                  and pv.dist2(g.vertices[j]).is_rational()
                  and pv.dist2(g.vertices[j]).c[0] >= Fraction(1, 4)]
        for b in others[:4]:
            copies = two_orbit_block(g, bp, legs[0], b)
            if copies is None:
                continue
            assert len(copies) == 6
            assert blocks_two_targets(g, bp, [legs[0], b], copies)
            assert cross_blocks(g, bp, [legs[0], b], copies)
            tested += 1
        if tested >= 12:
            break
    assert tested >= 6


def test_two_orbit_refuses_a_leg_off_the_classical_circle():
    from hn.mixed import two_orbit_block

    g = build_graph(build_Sa())
    pv = g.vertices[0]
    off = next(j for j in range(1, g.n) if j not in g.adj[0]
               and pv.dist2(g.vertices[j]).is_rational()
               and pv.dist2(g.vertices[j]).c[0] != __import__(
                   "fractions").Fraction(1, 3))
    try:
        two_orbit_block(g, 0, off, 1)
    except ValueError:
        return
    raise AssertionError("should refuse a leg that is not on d^2 = 1/3")


def test_centroids_sit_on_the_classical_circle():
    """A unit triangle's centroid is 1/sqrt(3) from each of its corners.

    Which is what manufactures the leg the two-orbit block needs. Points that
    close together are scarce in constructions built from unit steps, and on
    Sa the centroids take every tested pivot from having no leg to having one.
    """
    from fractions import Fraction

    from hn.mixed import unit_triangle_centroids

    g = build_graph(build_Sa())
    cents = unit_triangle_centroids(g, limit=60)
    assert cents
    third = Fraction(1, 3)
    from itertools import combinations

    for c in cents:
        close = [q for q in g.vertices
                 if c.dist2(q).is_rational() and c.dist2(q).c[0] == third]
        assert len(close) >= 3                 # its own triangle, at least
        # some triple among them is that triangle -- not necessarily the
        # first three, since a centroid can sit 1/sqrt(3) from more points
        # than the three it came from.
        assert any(u.dist2(v) == 1 and v.dist2(w) == 1 and u.dist2(w) == 1
                   for u, v, w in combinations(close, 3))


def test_centroids_give_every_pivot_a_leg():
    from fractions import Fraction

    from hn.mixed import unit_triangle_centroids

    g = build_graph(build_Sa())
    known = set(g.vertices)
    fresh = [c for c in unit_triangle_centroids(g) if c not in known]
    assert fresh
    g2 = build_graph(list(g.vertices) + fresh)
    third = Fraction(1, 3)
    has = 0
    for bp in range(min(g2.n, 40)):
        pv = g2.vertices[bp]
        if any(pv.dist2(g2.vertices[j]).is_rational()
               and pv.dist2(g2.vertices[j]).c[0] == third
               for j in range(g2.n) if j != bp):
            has += 1
    assert has == min(g2.n, 40)
