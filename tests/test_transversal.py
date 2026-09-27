"""The same-distance ceiling, and the measurement that corrected it."""
import pytest

from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD, rotation_joining
from hn.multispindle import rotation_powers
from hn.transversal import (
    cycle_lengths,
    diagnose,
    multiquadratic_orders,
    same_distance_ceiling,
    trapping_bound,
)


def test_multiquadratic_orders_are_the_divisors_of_24():
    assert multiquadratic_orders() == [1, 2, 3, 4, 6, 8, 12, 24]


def test_triangle_is_the_only_odd_cycle_available():
    odd = {q for n in multiquadratic_orders() for q in cycle_lengths(n) if q % 2}
    assert odd == {1, 3}


def test_bipartite_components_never_trap():
    for q in (2, 4, 6, 8, 12, 24):
        assert trapping_bound(q) == 0


def test_ceiling_is_two_and_odd_orders_lift_it():
    assert same_distance_ceiling(multiquadratic_orders()) == 2
    assert same_distance_ceiling(multiquadratic_orders() + [5]) == 3
    assert same_distance_ceiling(multiquadratic_orders() + [7]) == 4
    assert same_distance_ceiling(multiquadratic_orders() + [11]) == 5


@pytest.mark.parametrize("d2v", ["1/3", "5/9"])
def test_degree_two_among_distinct_images(d2v):
    """The claim the first version of this module got wrong.

    Counting degrees in the conflict graph suggested it grew without bound as
    copies were added.  That was multiplicity: among *distinct* points on one
    circle the degree is 2, because the unit circle about a point meets that
    circle twice.
    """
    from fractions import Fraction

    g = build_G()
    p, val = 0, Fraction(d2v)
    pv = g.vertices[p]
    targets = [j for j in range(g.n)
               if j != p and pv.dist2(g.vertices[j]).is_rational()
               and pv.dist2(g.vertices[j]).c[0] == val][:6]
    assert targets
    rot = rotation_joining(val, DEGREY_FIELD)
    images = [r.about(pv)(g.vertices[q])
              for r in rotation_powers(rot, 12) for q in targets]
    distinct = list({im: None for im in images})
    for a, x in enumerate(distinct):
        deg = sum(1 for b, y in enumerate(distinct) if a != b and x.is_unit_apart(y))
        assert deg <= 2


def test_diagnose_calls_the_stalled_core_hopeless():
    """Eleven same-distance targets against a ceiling of two."""
    g = build_G()
    pv = g.vertices[0]
    same = [j for j in range(1, g.n)
            if pv.dist2(g.vertices[j]) == pv.dist2(g.vertices[1])][:11]
    d = diagnose(g, 0, same)
    assert d["single_circle"] and d["ceiling"] == 2
    if len(same) > 2:
        assert d["hopeless"]


def test_magic_radius_is_the_famous_one_at_order_three():
    from hn.transversal import magic_radius

    assert abs(magic_radius(3) ** 2 - 1 / 3) < 1e-12


def test_capacity_hierarchy():
    """Computed, not read off the independence number, which over-counts."""
    from hn.transversal import largest_odd_divisor, trapping_bound, trapping_capacity

    assert largest_odd_divisor(24) == 3
    assert trapping_capacity(3) == 2          # every multiquadratic field
    assert trapping_capacity(24) == 2         # even at the largest order
    assert trapping_capacity(5) == 3
    assert [trapping_bound(q) for q in (3, 5, 7, 9, 11, 13, 15)] == [2, 3, 4, 4, 5, 5, 4]
    assert trapping_capacity(15) == 4         # its C_15 beats its C_5 and C_3


def test_capacity_peaks_at_five_and_eleven_is_out_of_reach():
    """The size the narrowing stalled at is past every capacity computed."""
    from hn.transversal import order_for_capacity, trapping_bound

    assert order_for_capacity(2) == 3
    assert order_for_capacity(4) == 7
    assert order_for_capacity(5) == 11
    assert order_for_capacity(6) == 0         # nothing up to length 15 holds 6
    assert max(trapping_bound(q) for q in range(3, 16, 2)) == 5


def test_niven_leaves_exactly_four_rational_radii():
    """cos t = 1 - 1/(2d^2) is rational iff d^2 is, and Niven allows four."""
    from fractions import Fraction

    from hn.transversal import niven_capacity, niven_order

    assert niven_order(1) == 6
    assert niven_order(Fraction(1, 2)) == 4
    assert niven_order(Fraction(1, 3)) == 3
    assert niven_order(Fraction(1, 4)) == 2
    for d2 in (2, 3, Fraction(5, 3), Fraction(1, 5), Fraction(2, 3), 7):
        assert niven_order(d2) == 0          # a path: nothing is trapped
        assert niven_capacity(d2) == 0


def test_only_the_classical_circle_traps_at_a_rational_radius():
    from fractions import Fraction

    from hn.transversal import niven_capacity

    assert niven_capacity(Fraction(1, 3)) == 2
    for d2 in (1, Fraction(1, 2), Fraction(1, 4)):
        assert niven_capacity(d2) == 0       # even orders are bipartite


# -- how large a core rotations can block ---------------------------------

def test_counting_certificate_matches_the_built_block():
    """alphas [2, 3] against 6 copies: the construction blocks by exactly one.

    The 1/3 leg's conflict graph is two triangles, so alpha = 2 = m/3, the
    theorem's floor; the other leg's is three disjoint edges, alpha = 3. Five
    is less than six, so counting alone rules out every escape, and the SAT
    test agrees.
    """
    from hn.graph import build_graph
    from hn.mixed import (blocks_targets, joint_core_configuration,
                          joint_core_copies)
    from hn.transversal import counting_blocks, leg_conflict_alphas

    E, pts, rho, sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    pair = [g.vertices.index(pts[2]), g.vertices.index(pts[3])]
    copies = joint_core_copies(E, rho, sigma)
    alphas = leg_conflict_alphas(g, bp, pair, copies)
    assert alphas == [2, 3] and sum(alphas) == 5 < len(copies) == 6
    assert counting_blocks(g, bp, pair, copies)
    assert blocks_targets(g, bp, pair, copies)


def test_every_leg_alpha_is_at_least_a_third_of_the_copies():
    """The floor the theorem rests on: images of one leg lie on one circle,
    a point of a circle is one apart from at most two points of it, so each
    conflict graph has maximum degree two and independence ratio >= 1/3.
    """
    from hn.graph import build_graph
    from hn.mixed import joint_core_configuration, joint_core_copies
    from hn.transversal import leg_conflict_alphas

    E, pts, rho, sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    legs = [g.vertices.index(pts[2]), g.vertices.index(pts[3])]
    copies = joint_core_copies(E, rho, sigma)
    for sub in (copies, copies[:3], copies[:5], copies[1:]):
        for a in leg_conflict_alphas(g, bp, legs, sub):
            assert 3 * a >= len(sub)


def test_counting_can_never_block_a_core_of_three():
    """r * m/3 <= sum alpha, so sum alpha < m forces r <= 2. Checked by
    adding a third leg to the built configuration: the sum reaches the copy
    count and the counting certificate dies, whatever the copies."""
    from hn.graph import build_graph
    from hn.mixed import joint_core_configuration, joint_core_copies
    from hn.transversal import counting_blocks, leg_conflict_alphas

    E, pts, rho, sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    trio = [g.vertices.index(pts[1]), g.vertices.index(pts[2]),
            g.vertices.index(pts[3])]
    copies = joint_core_copies(E, rho, sigma)
    for sub in (copies, copies[:4], copies[:5]):
        alphas = leg_conflict_alphas(g, bp, trio, sub)
        assert sum(alphas) >= len(sub), "counting must fail at three legs"
        assert not counting_blocks(g, bp, trio, sub)


def test_three_legs_block_and_four_do_not():
    """Misalignment reaches three and stops; the counting slack says why.

    Each leg's conflict graph has maximum degree two, so the independent sets
    available total r*m/3. At three legs that is exactly the copy count -- the
    knife edge where overlap still decides it -- and at four it is a third
    more room than there are copies to place.
    """
    from pysat.solvers import Solver

    from hn.transversal import MAX_BLOCKABLE_CORE

    def blocks(N, ts, ss):
        r = len(ts)

        def x(k, e, L):
            return 1 + ((k * 2 + e) * r + L)

        cls = [[x(k, e, L) for L in range(r)]
               for k in range(N) for e in (0, 1)]
        for L in range(r):
            t, s = ts[L], ss[L]
            for k in range(N):
                for e in (0, 1):
                    for d in (t, -t):
                        w = (k - d) % N
                        if (w, e) != (k, e):
                            cls.append([-x(k, e, L), -x(w, e, L)])
                    for d in (t - s, -t - s):
                        cls.append([-x(k, 0, L), -x((k - d) % N, 1, L)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            return not sv.solve()

    assert MAX_BLOCKABLE_CORE == 3
    # the smallest known three-leg block: all legs at d^2 = 1/3
    assert blocks(3, (1, 1, 1), (2, 2, 0))
    # rotations only -- every shift zero -- never blocks three
    assert not blocks(3, (1, 1, 1), (0, 0, 0))
    # and four legs do not block, on the same N and angles, whatever the shifts
    for s3 in range(3):
        for s4 in range(3):
            assert not blocks(3, (1, 1, 1, 1), (2, 2, s3, s4))
