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
    assert odd_cycle_reach(15) == 8


def test_arithmetic_is_consistent_with_the_complex_embedding():
    R = CycloRing(15)
    us = R.unit_steps()
    a, b = us[1], us[4]
    prod = R.mul(a, b)
    assert abs(R.to_complex(prod) - R.to_complex(a) * R.to_complex(b)) < 1e-9
    diff = R.sub(a, b)
    assert abs(abs(R.to_complex(diff)) ** 2
               - R.to_complex(R.norm2(diff)).real) < 1e-9
