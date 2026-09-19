"""The field is the foundation: if arithmetic here is wrong, every edge is."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fractions import Fraction

import pytest

from hn.field import Field, QSQRT3_11


def test_generators_must_be_squarefree_and_coprime():
    with pytest.raises(ValueError):
        Field((4, 11))       # 4 = 2^2
    with pytest.raises(ValueError):
        Field((6, 15))       # share the factor 3


def test_radicals_square_to_their_integers():
    F = QSQRT3_11
    for d in (3, 11, 33):
        assert F.sqrt(d) * F.sqrt(d) == d


def test_product_of_radicals():
    F = QSQRT3_11
    assert F.sqrt(3) * F.sqrt(11) == F.sqrt(33)


def test_sqrt_extracts_square_factors():
    F = QSQRT3_11
    assert F.sqrt(12) == F.sqrt(3) * 2
    assert F.sqrt(99) == F.sqrt(11) * 3
    with pytest.raises(ValueError):
        F.sqrt(5)            # would need a larger field


def test_inverse_uses_every_galois_conjugate():
    F = QSQRT3_11
    x = F.rational(Fraction(5, 6)) + F.sqrt(11) * F.rational(Fraction(1, 6))
    assert x * x.inverse() == 1
    y = F.rational(3) + F.sqrt(3) * 2 - F.sqrt(11) + F.sqrt(33) * Fraction(7, 5)
    assert y * y.inverse() == 1
    assert (y / y) == 1


def test_zero_has_no_inverse():
    with pytest.raises(ZeroDivisionError):
        QSQRT3_11.zero().inverse()


def test_powers_and_negative_powers():
    F = QSQRT3_11
    x = F.rational(Fraction(5, 6)) + F.sqrt(11) * F.rational(Fraction(1, 6))
    assert x ** 3 == x * x * x
    assert (x ** -2) * (x ** 2) == 1


def test_equality_is_exact_not_approximate():
    F = QSQRT3_11
    # 1.7320508... vs a good rational approximation: floats would call these equal
    approx = F.rational(Fraction(1732050807, 1000000000))
    assert F.sqrt(3) != approx
    assert abs(float(F.sqrt(3)) - float(approx)) < 1e-9


def test_hash_matches_equality():
    F = QSQRT3_11
    a = F.sqrt(3) * 2 + F.rational(1)
    b = F.rational(1) + F.sqrt(3) * 2
    assert a == b and hash(a) == hash(b)
