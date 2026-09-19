"""The general spindle lemma, and the arithmetic that bounds it."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fractions import Fraction
from math import gcd

import pytest

from hn.coloring import is_k_colorable
from hn.field import Field, QSQRT3_11
from hn.geometry import DEGREY_FIELD, ROT60, Point, eisenstein, origin, rotation_joining
from hn.graph import build_graph
from hn.multispindle import (MULTIQUADRATIC_ORDERS, available_rotation_orders,
                             blocks_all_assignments, conflict_graphs,
                             independent_transversal, multi_spindle_union,
                             rotation_powers, spindle_catalogue,
                             squared_distance_for_step)

FIELDS = [QSQRT3_11, DEGREY_FIELD, Field((2, 3, 5, 7, 11)), Field((2, 3))]


# -- the lemma reproduces what it should ----------------------------------

def test_two_copy_case_is_the_classic_spindle():
    rhombus = build_graph([origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)])
    rots = rotation_powers(rotation_joining(3), 2)
    assert blocks_all_assignments(rhombus, 0, [3], rots)
    spun = multi_spindle_union(rhombus, 0, rots)
    assert (spun.n, spun.m) == (7, 11)
    assert is_k_colorable(spun, 3)[0] is False


def unit_triangle():
    r = QSQRT3_11.sqrt(3).inverse()
    q = Point(r, QSQRT3_11.zero())
    rho = ROT60 ** 2
    return build_graph([origin(), q, rho(q), rho(rho(q))])


def test_three_copies_block_two_targets_but_not_three():
    g = unit_triangle()
    rots = rotation_powers(ROT60 ** 2, 3)
    assert blocks_all_assignments(g, 0, [1, 2], rots)
    # three copies cannot block three targets: the bijection is a transversal
    assert not blocks_all_assignments(g, 0, [1, 2, 3], rots)
    assert independent_transversal(3, [1, 2, 3], conflict_graphs(g, 0, [1, 2, 3], rots))


# -- the arithmetic ceiling ------------------------------------------------

def test_realisable_orders_are_exactly_the_divisors_of_24():
    """A rotation of order n puts zeta_n in F(i); when F is multiquadratic
    that forces (Z/n)* to have exponent 2, i.e. n | 24."""
    # n = 1 is degenerate (x*x % 1 == 0), so the criterion is applied from 2 up
    predicted = [1] + [n for n in range(2, 200)
                       if all((x * x) % n == 1 for x in range(1, n + 1) if gcd(x, n) == 1)]
    assert predicted == [1, 2, 3, 4, 6, 8, 12, 24]
    assert set(MULTIQUADRATIC_ORDERS) == set(predicted)


def test_three_is_the_only_odd_order():
    assert [n for n in MULTIQUADRATIC_ORDERS if n % 2 and n > 1] == [3]


@pytest.mark.parametrize("field", FIELDS)
def test_orders_need_their_radicals(field):
    orders = available_rotation_orders(field)
    assert 3 in orders if 3 in field.gens else True
    # 8 and 24 need sqrt2; 3, 6 and 12 need sqrt3
    assert (8 in orders) == (2 in field.gens)
    assert (24 in orders) == (2 in field.gens and 3 in field.gens)
    assert (12 in orders) == (3 in field.gens)


@pytest.mark.parametrize("field", FIELDS)
def test_catalogue_reaches_no_further_than_two_targets(field):
    """The reach of the CATALOGUE, which is not the same as a theorem.

    spindle_catalogue only searches cyclic families {rho^0..rho^(n-1)}.  Over
    those it finds nothing with three targets.  Arbitrary finite families --
    a torsion rotation mixed with an infinite-order one, say -- are allowed by
    the lemma and are not covered here, so this pins current reach, not
    impossibility.
    """
    cat = spindle_catalogue(field, max_targets=3)
    assert cat, "field supports no spindles at all"
    assert all(len(e["steps"]) <= 2 for e in cat)


def test_catalogue_contains_the_known_triple_spindle():
    cat = spindle_catalogue(QSQRT3_11, max_targets=2)
    hit = [e for e in cat
           if e["n"] == 3 and len(e["steps"]) == 2 and all(d == Fraction(1, 3) for d in e["d2"])]
    assert hit, "the classic d^2 = 1/3 three-copy spindle is missing"


def test_bipartite_conflict_graphs_can_still_block_two_targets():
    """Order 4 has no odd cycle, yet steps (1, 2) block two targets.

    So "blocking needs a non-bipartite conflict graph" is false once the two
    graphs differ; chromatic number governs only the all-equal case.
    """
    cat = spindle_catalogue(QSQRT3_11, max_targets=2)
    mixed = [e for e in cat if len(e["steps"]) == 2 and e["d2"][0] != e["d2"][1]]
    assert len(mixed) >= 5
    pairs = {(e["n"], str(e["d2"][0]), str(e["d2"][1])) for e in mixed}
    assert (4, "1/2", "1/4") in pairs      # order 4: no odd cycle anywhere
    assert (6, "1/3", "1/4") in pairs


def test_adding_sqrt2_widens_the_catalogue():
    small = spindle_catalogue(QSQRT3_11, max_targets=2)
    big = spindle_catalogue(Field((2, 3, 11)), max_targets=2)
    assert len(big) > len(small)
    assert 8 in available_rotation_orders(Field((2, 3, 11)))
    assert 8 not in available_rotation_orders(QSQRT3_11)


def test_targets_at_distance_one_are_rejected():
    """Such a target is adjacent to the pivot, so it can never be the one
    sharing its colour, and the step is useless."""
    assert squared_distance_for_step(6, 1, QSQRT3_11) is None      # 60 deg -> d^2 = 1
    assert squared_distance_for_step(6, 2, QSQRT3_11) == Fraction(1, 3)
    assert squared_distance_for_step(4, 2, QSQRT3_11) == Fraction(1, 4)


def test_catalogued_distances_are_exact_field_elements():
    """Each entry carries its own rotations, since a catalogued d^2 may be
    irrational (2 + sqrt(3) at order 12) and rotation_joining takes rationals."""
    for e in spindle_catalogue(QSQRT3_11, max_targets=2):
        assert len(e["rotations"]) == e["n"]
        for d2 in e["d2"]:
            assert d2.field == QSQRT3_11
            assert not d2.is_zero()
        for step, d2 in zip(e["steps"], e["d2"]):
            rot = e["rotations"][step]
            # a point at squared distance d2 really does move by exactly 1
            assert (rot.cos - 1) * (rot.cos - 1) + rot.sin * rot.sin == (d2 * 2).inverse() * 2


# -- how much forcing machinery a field actually offers ---------------------

def test_spindle_spectrum_grows_with_the_field():
    """Counting spindle-able distances measures the search space a field gives.

    The literature works in Q(sqrt3, sqrt11); de Grey's rotations force the
    larger field; adding sqrt2 -- which nothing in the sources read here does --
    widens it further again.
    """
    from hn.multispindle import spindle_spectrum

    small = spindle_spectrum(QSQRT3_11)
    degrey = spindle_spectrum(DEGREY_FIELD)
    with_root2 = spindle_spectrum(Field((2, 3, 5, 7, 11)))
    assert len(small) < len(degrey) < len(with_root2)
    assert set(small) <= set(degrey) <= set(with_root2)


def test_spectrum_entries_really_have_spindle_rotations():
    from hn.multispindle import spindle_spectrum

    for d2 in spindle_spectrum(QSQRT3_11)[:20]:
        rot = rotation_joining(d2, QSQRT3_11)
        assert rot.cos * rot.cos + rot.sin * rot.sin == 1


def test_a_distance_outside_the_spectrum_is_rejected():
    from hn.multispindle import spindle_spectrum

    spectrum = set(spindle_spectrum(QSQRT3_11))
    assert Fraction(25, 9) not in spectrum        # needs sqrt91 = sqrt7*sqrt13
    with pytest.raises(ValueError):
        rotation_joining(Fraction(25, 9), QSQRT3_11)
    # and it does exist once the field carries those radicals
    rotation_joining(Fraction(25, 9), Field((7, 13)))
