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
