"""Z[zeta_n]: the unit steps, the triangles, and what the field buys."""
import pytest

from hn.cyclotomic import CycloRing, cyclotomic_poly, odd_cycle_reach, rotation_orders


def test_cyclotomic_polynomials():
    assert cyclotomic_poly(1) == (-1, 1)
    assert cyclotomic_poly(3) == (1, 1, 1)
    assert cyclotomic_poly(4) == (1, 0, 1)
    assert cyclotomic_poly(5) == (1, 1, 1, 1, 1)
    assert len(cyclotomic_poly(15)) - 1 == 8          # phi(15)


@pytest.mark.parametrize("n,steps", [(3, 6), (4, 4), (5, 10), (15, 30)])
def test_unit_steps_are_the_roots_of_unity(n, steps):
    R = CycloRing(n)
    us = R.unit_steps()
    assert len(us) == steps
    for z in us:
        assert R.norm2(z) == R.one()


def test_eisenstein_is_the_n_equals_three_case():
    """Six unit steps and unit triangles: the world every search here used."""
    R = CycloRing(3)
    us = R.unit_steps()
    assert len(us) == 6
    tri = [1 for i, a in enumerate(us) for j, b in enumerate(us) if j > i
           for c in us if R.add(R.add(a, b), c) == R.zero()]
    assert tri


def test_zeta_five_alone_has_no_unit_triangle():
    """Three 10th roots of unity never sum to zero, so the clique number drops."""
    R = CycloRing(5)
    us = R.unit_steps()
    tri = [1 for i, a in enumerate(us) for j, b in enumerate(us) if j > i
           for c in us if R.add(R.add(a, b), c) == R.zero()]
    assert not tri


def test_zeta_fifteen_keeps_the_triangle_and_lifts_the_ceiling():
    R = CycloRing(15)
    us = R.unit_steps()
    tri = [1 for i, a in enumerate(us) for j, b in enumerate(us) if j > i
           for c in us if R.add(R.add(a, b), c) == R.zero()]
    assert tri                                        # omega is still in there
    assert 5 in rotation_orders(15) and 15 in rotation_orders(15)
    assert odd_cycle_reach(3) == 2
    assert odd_cycle_reach(15) == 4


def test_arithmetic_is_consistent_with_the_complex_embedding():
    R = CycloRing(15)
    us = R.unit_steps()
    a, b = us[1], us[4]
    prod = R.mul(a, b)
    assert abs(R.to_complex(prod) - R.to_complex(a) * R.to_complex(b)) < 1e-9
    diff = R.sub(a, b)
    assert abs(abs(R.to_complex(diff)) ** 2
               - R.to_complex(R.norm2(diff)).real) < 1e-9


def test_gauss_sum_gives_the_spindle_radical():
    from hn.cyclotomic import CycloField, moser_rotation

    F = CycloField(33)
    r = F.sqrt_disc(11)
    assert F.mul(r, r) == F.rational(-11)
    rho = moser_rotation(F)
    assert F.norm2(rho) == F.one()          # a rotation, so modulus one


def test_moser_spindle_lives_in_a_cyclotomic_field():
    """Four-chromatic, over a field with odd-order rotations available.

    The spindle's rotation is multiplication by (5 + sqrt(-11))/6, and
    sqrt(-11) is the Gauss sum over zeta_11, so the whole construction fits in
    Q(zeta_33) -- which also carries rotations of order 11 and 33, capacities
    5 and 5 against the multiquadratic 2.
    """
    from hn.coloring import is_k_colorable
    from hn.cyclograph import build_cyclo_graph, moser_spindle_cyclotomic
    from hn.cyclotomic import CycloField
    from hn.transversal import trapping_capacity

    F = CycloField(33)
    g = build_cyclo_graph(moser_spindle_cyclotomic(F))
    assert (g.n, g.m) == (7, 11)
    assert not is_k_colorable(g, 3)[0]
    assert is_k_colorable(g, 4)[0]
    assert trapping_capacity(11) == 5 and trapping_capacity(33) == 5


def test_unit_polygon_realises_the_magic_circle():
    """The magic circle of order n is the circumcircle of a unit n-gon.

    And the n-gon has a closed form with no radicals: v_k = z^k / (z - 1),
    whose consecutive vertices are |z - 1| / |z - 1| = 1 apart. For odd n it
    is an odd unit-distance cycle, which is the only structure that traps.
    """
    import math

    from hn.cyclotomic import CycloField, unit_polygon

    F = CycloField(33)
    for order in (3, 11, 33):
        vs = unit_polygon(F, order)
        assert len(vs) == order
        for k in range(order):
            assert F.norm2(F.sub(vs[(k + 1) % order], vs[k])) == F.one()
        r2 = F.to_complex(F.norm2(vs[0])).real
        assert abs(r2 - 1 / (4 * math.sin(math.pi / order) ** 2)) < 1e-9


def test_the_classical_magic_circle_is_the_order_three_case():
    from fractions import Fraction

    from hn.cyclotomic import CycloField, unit_polygon

    F = CycloField(33)
    v = unit_polygon(F, 3)[0]
    assert F.norm2(v) == F.rational(Fraction(1, 3))
