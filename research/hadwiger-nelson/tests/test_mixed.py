"""Mixed-distance machinery: denesting, conflict rotations, escape counts."""
from fractions import Fraction

from hn.degrey import build_Sa
from hn.field import Field
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.certify import exact_edges
from hn.coloring import is_k_colorable
from hn.mixed import (_sqrt_in_field, blocks_two_targets,
                      compose_rotations, conflict_rotation_set,
                      count_cross_transversals, joint_core_configuration,
                      joint_core_copies, joint_core_union)
from hn.spindle import SeparationTest


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


def test_two_orbit_block_runs_to_an_uncolourable_union():
    """The whole chain on a real object: forced pair, six copies, no colouring.

    What it shows: the machinery runs end to end and the union really has no
    colouring. 91 vertices go in 3-colourable; the union of the six copies
    comes out at 409 vertices and 1062 edges with none.

    What it does *not* show: that the block reaches anything the classical
    argument cannot. The sqrt(3) leg of this pair is forced on its own here --
    a core of 1, since two points at distance sqrt(3) in the triangular
    lattice take the same colour in every 3-colouring -- so the pair is forced
    trivially and two copies of the ordinary spindle would close it too. The
    test asserts that, so the limitation cannot quietly disappear.

    A genuine demonstration needs both legs separable alone. That object has
    not turned up yet, at any k.
    """
    from fractions import Fraction

    from hn.coloring import is_k_colorable
    from hn.generate import hex_ball
    from hn.mixed import (blocks_two_targets, two_orbit_block,
                          unit_triangle_centroids)
    from hn.multispindle import cross_blocks
    from hn.spindle import SeparationTest

    pts = list(hex_ball(3))
    g0 = build_graph(pts)
    cents = [c for c in unit_triangle_centroids(g0) if c not in set(g0.vertices)]
    g = build_graph(pts + cents)
    assert is_k_colorable(g, 3)[0]

    third = Fraction(1, 3)
    found = None
    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v])):
        p = g.vertices[bp]
        legs, others = [], []
        for j in range(g.n):
            if j == bp or j in g.adj[bp]:
                continue
            d2 = p.dist2(g.vertices[j])
            if not d2.is_rational():
                continue
            (legs if d2.c[0] == third else others).append(j)
        if not legs or not others:
            continue
        st = SeparationTest(g, 3, bp, legs + others)
        try:
            for a in legs:
                for b in others:
                    if not st.run(subset=[a, b])[0]:
                        found = (bp, a, b)
                        break
                if found:
                    break
        finally:
            st.close()
        if found:
            break
    assert found, "no forced pair with a leg on the classical circle"

    bp, a, b = found
    st = SeparationTest(g, 3, bp, [a, b])
    try:
        # pinned deliberately: on this instance the pair is forced only
        # because one leg already is, and that is the whole limitation
        assert not st.run(subset=[b])[0]
    finally:
        st.close()
    copies = two_orbit_block(g, bp, a, b)
    assert copies is not None and len(copies) == 6
    assert blocks_two_targets(g, bp, [a, b], copies)
    assert cross_blocks(g, bp, [a, b], copies)

    p = g.vertices[bp]
    union = {q: None for q in g.vertices}
    for r in copies:
        f = r.about(p)
        for q in g.vertices:
            union[f(q)] = None
    u = build_graph(list(union))
    assert not is_k_colorable(u, 3)[0]


def test_odd_orbit_block_generalises_the_six_copy_one():
    """Two orbits of any odd-order rotation close a forced pair, with 2n copies.

    Nothing in the six-copy argument needed the triangle, only that the cycle
    be odd. An orbit's images of the leg form C_order(step), whose independent
    sets hold at most (order-1)/2 elements, so at least (order+1)/2 copies per
    orbit must choose the other target -- and two subsets of {0..order-1} of
    that size total at least order+1, so they intersect.

    At order 3 it must agree with two_orbit_block exactly, which is what this
    pins. Orders 5 and up need a rotation of that order, and a multiquadratic
    field has none: n | 24 leaves 1 and 3 as its only odd orders, so the
    generalisation buys nothing here and everything in a cyclotomic field.
    """
    from fractions import Fraction

    from hn.mixed import blocks_two_targets, odd_orbit_block, two_orbit_block
    from hn.multispindle import cross_blocks

    g = build_graph(build_Sa())
    third = Fraction(1, 3)
    checked = 0
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
        for b in others[:3]:
            old = two_orbit_block(g, bp, legs[0], b)
            new = odd_orbit_block(g, bp, legs[0], b, order=3)
            if old is None or new is None:
                continue
            assert len(new) == 6
            assert blocks_two_targets(g, bp, [legs[0], b], new)
            assert cross_blocks(g, bp, [legs[0], b], new)
            checked += 1
        if checked >= 6:
            break
    assert checked >= 3


def test_odd_orbit_block_refuses_even_or_non_coprime_orders():
    from hn.mixed import odd_orbit_block

    g = build_graph(build_Sa())
    for order, step in ((4, 1), (2, 1), (9, 3)):
        try:
            odd_orbit_block(g, 0, 1, 2, order=order, step=step)
        except ValueError:
            continue
        raise AssertionError(f"should refuse order={order}, step={step}")


def test_a_genuine_core_of_two_exists_on_four_points():
    """The object the whole search is for, built by hand from K_4.

    A pair {a,b} is forced exactly when the graph plus the two edges (p,a),
    (p,b) has no k-colouring -- so the smallest instance is the smallest
    (k+1)-critical graph with two edges removable at one vertex, which at
    k = 3 is K_4.

    Take a unit triangle x, y, z and a point p at distance 1 from x and
    1/sqrt(3) from y. The triangle consumes all three colours and p differs
    from x, so p must share with y or with z: forced. Neither leg is forced
    alone, since p may take either. And y sits on the classical circle.

    Its geometry is rigid: |p - x| = 1 and |p - y|^2 = 1/3 leave two positions
    for p, and both give |p - z|^2 = (7 +- sqrt(33))/6. One is 0.209, too
    short for a unit chord; the other is 2.124, long enough, but the rotation
    it needs has sin^2 = (-66 + 30 sqrt(33))/256, whose conjugate is negative
    -- so no real multiquadratic field names it. The pair is genuine; the
    field is what stops the block being applied to it.
    """
    from fractions import Fraction

    from hn.field import Field
    from hn.geometry import Point
    from hn.spindle import SeparationTest

    fd = Field((3, 11))
    x = Point(fd.rational(0), fd.rational(0))
    y = Point(fd.rational(1), fd.rational(0))
    z = Point(fd.rational(Fraction(1, 2)),
              fd.sqrt(3) * fd.rational(Fraction(1, 2)))
    p = Point(fd.rational(Fraction(5, 6)),
              fd.sqrt(11) * fd.rational(Fraction(-1, 6)))

    assert x.dist2(y) == 1 and y.dist2(z) == 1 and x.dist2(z) == 1
    assert p.dist2(x) == 1
    assert p.dist2(y) == fd.rational(Fraction(1, 3))
    assert p.dist2(z) == fd.rational(Fraction(7, 6)) + fd.sqrt(33) * fd.rational(
        Fraction(1, 6))

    g = build_graph([p, x, y, z])
    idx = {q: i for i, q in enumerate(g.vertices)}
    bp, a, b = idx[p], idx[y], idx[z]
    assert set(g.adj[bp]) == {idx[x]}            # p meets only x

    st = SeparationTest(g, 3, bp, [a, b])
    try:
        assert st.run(subset=[a])[0]             # y alone: separable
        assert st.run(subset=[b])[0]             # z alone: separable
        assert not st.run(subset=[a, b])[0]      # the pair: forced
    finally:
        st.close()

    from hn.mixed import two_orbit_block
    assert two_orbit_block(g, bp, a, b) is None  # the field cannot name sigma


# -- the genuinely joint core of two --------------------------------------
#
# Every forced pair found by search was forced one leg at a time; the
# 409-vertex certificate blocks a pair whose 1/3 leg is already forced alone,
# so the block was valid and unnecessary. This configuration is the case the
# classical argument cannot reach: neither leg forced by itself, the pair
# forced only jointly, and six copies closing it. These tests pin the four
# points, the two rotations, the joint forcing, and the union's chromatic
# number, because every one of them is load-bearing for the claim.

def _genuine_pair():
    E, pts, rho, sigma = joint_core_configuration()
    return E, pts, rho, sigma, compose_rotations


def test_genuine_pair_distances_are_exact():
    E, (p, x, y, z), _rho, _sigma, _c = _genuine_pair()
    one = E.rational(1)
    assert p.dist2(x) == one                      # p is adjacent to x
    assert p.dist2(y) == E.rational(Fraction(1, 3))
    # (7 + sqrt 33)/6, the leg no multiquadratic rotation can close
    base = E.base
    assert p.dist2(z) == E.embed(base.rational(Fraction(7, 6))
                                 + base.sqrt(33) * base.rational(Fraction(1, 6)))
    assert x.dist2(y) == one and y.dist2(z) == one and z.dist2(x) == one


def test_genuine_pair_core_is_two_and_joint():
    """Neither leg alone, both together. This is the whole point."""
    E, (p, x, y, z), _rho, _sigma, _c = _genuine_pair()
    g = build_graph([p, x, y, z])
    bp = g.vertices.index(p)
    ia, ib = g.vertices.index(y), g.vertices.index(z)
    st = SeparationTest(g, 3, bp, [ia, ib])
    try:
        assert st.run(subset=[ia])[0], "y alone must NOT be forced"
        assert st.run(subset=[ib])[0], "z alone must NOT be forced"
        assert not st.run(subset=[ia, ib])[0], "the pair must be forced"
    finally:
        st.close()


def test_genuine_pair_rotations_close_their_circles():
    E, (p, _x, y, z), rho, sigma, _c = _genuine_pair()
    one = E.rational(1)
    assert sigma.cos * sigma.cos + sigma.sin * sigma.sin == one
    assert y.dist2(rho.about(p)(y)) == one        # 120 degrees on the 1/3 circle
    assert z.dist2(sigma.about(p)(z)) == one      # sqrt(v) on the other


def test_genuine_pair_union_is_not_three_colourable():
    E, (p, x, y, z), rho, sigma, compose = _genuine_pair()
    copies = joint_core_copies(E, rho, sigma)
    g4 = build_graph([p, x, y, z])
    bp = g4.vertices.index(p)
    assert blocks_two_targets(g4, bp, [g4.vertices.index(y),
                                       g4.vertices.index(z)], copies)
    pts = joint_core_union()
    g = build_graph(pts)
    assert (g.n, g.m) == (19, 33)
    assert set(exact_edges(pts)) == set(g.edges())
    assert not is_k_colorable(g, 3)[0]
    assert is_k_colorable(g, 4)[0]


def test_genuine_pair_union_is_vertex_critical():
    """No vertex is spare: the block uses every copy it builds."""
    pts = joint_core_union()
    for i in range(len(pts)):
        rest = build_graph([q for j, q in enumerate(pts) if j != i])
        assert is_k_colorable(rest, 3)[0], f"vertex {i} was not needed"
