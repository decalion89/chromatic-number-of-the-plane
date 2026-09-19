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
