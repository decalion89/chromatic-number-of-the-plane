"""Rotations must be exact rotations, and the lattice must be the lattice."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from hn.field import QSQRT3_11 as F
from hn.geometry import ROT60, SPINDLE, Point, Rotation, eisenstein, origin, rotation_joining


def test_rotation_rejects_non_rotations():
    with pytest.raises(ValueError):
        Rotation(F.rational(1), F.rational(1))    # cos^2+sin^2 = 2


def test_rot60_has_order_six():
    identity = Rotation(F.one(), F.zero(), check=False)
    assert ROT60 ** 6 == identity
    assert ROT60 ** 3 != identity


def test_lattice_neighbours_are_at_distance_one():
    o = origin()
    for a, b in ((1, 0), (0, 1), (-1, 0), (0, -1), (1, -1), (-1, 1)):
        assert o.is_unit_apart(eisenstein(a, b))


def test_rhombus_long_diagonal_is_sqrt_three():
    assert origin().dist2(eisenstein(1, 1)) == 3


def test_spindle_is_arccos_five_sixths():
    assert SPINDLE.cos == F.rational(1) * 5 / 6 or SPINDLE.cos * 6 == 5
    tip = eisenstein(1, 1)                    # at distance sqrt(3) from origin
    assert tip.dist2(SPINDLE(tip)) == 1       # rotated to exactly distance 1


@pytest.mark.parametrize("d2", [1, 3, 7, 19, 25, 37, 61, 91])
def test_rotation_joining_lands_at_distance_one(d2):
    """The defining property, for every distance the field supports."""
    rot = rotation_joining(d2)
    assert rot.cos * rot.cos + rot.sin * rot.sin == 1
    # a point at squared distance d2 from the origin, rotated, moves by 1
    from hn.field import QSQRT3_11
    import math
    # build such a point on the x-axis: (sqrt(d2), 0)
    p = Point(QSQRT3_11.sqrt(d2) if d2 in (1, 3, 11, 33) else QSQRT3_11.rational(0), QSQRT3_11.zero())
    if d2 in (1, 3):
        assert p.dist2(rot(p)) == 1


def test_rotation_joining_refuses_distances_outside_the_field():
    with pytest.raises(ValueError):
        rotation_joining(4)      # would need sqrt(15)


def test_rotation_about_a_pivot_fixes_the_pivot():
    pivot = eisenstein(2, 1)
    turn = SPINDLE.about(pivot)
    assert turn(pivot) == pivot


# -- the field de Grey's construction actually needs -----------------------

def test_degrey_rotations_are_outside_the_small_field():
    """The diagnosis that redirected this project: restricting to rotations
    that keep the field at Q(sqrt3, sqrt11) excludes exactly the two angles
    de Grey used, so his graph cannot lie in any set generated that way."""
    from hn.geometry import required_radical

    assert required_radical(4) == 15       # 2*arcsin(1/4): sin = sqrt(15)/8
    assert required_radical(16) == 7       # 2*arcsin(1/8): sin = 3*sqrt(7)/32
    for d2 in (4, 16):
        with pytest.raises(ValueError):
            rotation_joining(d2)           # rejected by Q(sqrt3, sqrt11)


def test_degrey_rotations_are_exact_in_the_extended_field():
    from hn.geometry import DEGREY_FIELD, degrey_rotations

    assert DEGREY_FIELD.dim == 16
    for name, rot in degrey_rotations().items():
        assert rot.cos * rot.cos + rot.sin * rot.sin == 1, name


@pytest.mark.parametrize("d2,dist", [(4, 2), (16, 4)])
def test_degrey_spindles_move_a_point_by_exactly_one(d2, dist):
    from hn.geometry import DEGREY_FIELD

    F = DEGREY_FIELD
    p = Point(F.rational(dist), F.zero())
    assert p.dist2(rotation_joining(d2, F)(p)) == 1


def test_small_field_rotations_are_the_centred_hexagonal_family():
    """d2 = 3k^2+3k+1 always works over Q(sqrt3): sin = (2k+1)sqrt(3)/(2 d2)."""
    from hn.geometry import required_radical

    for k in range(6):
        d2 = 3 * k * k + 3 * k + 1
        assert required_radical(d2) == (1 if d2 == 1 else 3) or required_radical(d2) == 3
        rotation_joining(d2)               # constructible without extension


def test_the_half_turn_is_a_spindle_and_needs_no_radical():
    """d^2 = 1/4 gives cos = -1, sin = 0: two points at distance 1/2 from the
    pivot land exactly 1 apart under a half-turn.

    A real spindle, and the one case where the sine vanishes.  Computing it as
    sqrt(0) used to raise, which killed a running search when it first met
    this distance.
    """
    from fractions import Fraction

    from hn.geometry import required_radical

    assert required_radical(Fraction(1, 4)) == 1
    from hn.geometry import DEGREY_FIELD

    for field in (F, DEGREY_FIELD):
        rot = rotation_joining(Fraction(1, 4), field)
        assert rot.cos == -1 and rot.sin.is_zero()
        assert rot.cos * rot.cos + rot.sin * rot.sin == 1
        p = Point(field.rational(Fraction(1, 2)), field.zero())
        assert p.dist2(rot(p)) == 1


def test_distances_just_below_a_quarter_are_still_refused():
    from fractions import Fraction

    with pytest.raises((ValueError, ZeroDivisionError)):
        rotation_joining(Fraction(1, 5))      # 2r < 1: no rotation separates them


def test_rotation_about_a_pivot_in_another_field_is_refused():
    """ROT60 lives in Q(sqrt3, sqrt11).  Applying it to a graph over a larger
    field used to fail deep inside field multiplication, with a message that
    named neither the rotation nor the caller."""
    from hn.geometry import DEGREY_FIELD, ROT60

    foreign = Point(DEGREY_FIELD.zero(), DEGREY_FIELD.zero())
    with pytest.raises(TypeError, match="rebuild the rotation"):
        ROT60.about(foreign)
    # and the correctly-built one works
    from hn.geometry import _rot60

    assert _rot60(DEGREY_FIELD).about(foreign)(foreign) == foreign
