"""Q(zeta_15, sqrt(-7), sqrt(-11)): de Grey's rotations and zeta_15 together."""
from fractions import Fraction

from hn.cyclotomic import CycloField
from hn.quadext import QuadExtField, degrey_field, degrey_rotations


def test_quadratic_extension_arithmetic():
    F = QuadExtField(CycloField(15), -11)
    assert F.degree == 16
    r = F.radical()
    assert F.mul(r, r) == F.rational(-11)
    assert F.conj(r) == F.neg(r)                 # purely imaginary
    assert F.norm2(F.moser_rotation()) == F.one()


def test_extensions_nest():
    L = degrey_field()
    assert L.degree == 32
    assert L.base.degree == 16 and L.base.base.degree == 8


def test_degrey_rotations_are_rotations_here():
    """Read as complex numbers they need sqrt(-15), sqrt(-7), sqrt(-11).

    As (cos, sin) pairs they look like they need the real multiquadratic field
    Q(sqrt3, sqrt5, sqrt7, sqrt11). They do not: that is only where the
    coordinates were written down.
    """
    L = degrey_field()
    rots = degrey_rotations(L)
    assert set(rots) == {"spindle_2", "spindle_4", "quarter", "moser"}
    for z in rots.values():
        assert L.norm2(z) == L.one()


def test_zeta_fifteen_is_in_degreys_field():
    """Which is the point: capacity 4 and his construction in one place."""
    from hn.transversal import trapping_bound

    L = degrey_field()
    z15 = L.embed(L.base.embed(L.base.base.zeta(1)))
    assert L.norm2(z15) == L.one()
    power = L.one()
    for _ in range(15):
        power = L.mul(power, z15)
    assert power == L.one()
    assert trapping_bound(15) == 4 and trapping_bound(3) == 2
