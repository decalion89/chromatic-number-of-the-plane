"""The unit-distance graph on the Moser field Q(sqrt-3, sqrt-11) is 4-colourable (hn/adelic.py).

2 does not split from Q(sqrt33) to Q(sqrt-3, sqrt-11), so every unit vector of the field is a
norm-one unit of Z_2[omega], and phi(alpha + beta omega) = frac_2((alpha + 2 beta)/4) is 1/4, 1/2 or
3/4 on all of them.  With the Moser spindle inside the field, its chromatic number is exactly 4.
"""
import itertools
import json
from fractions import Fraction as Fr

from hn.adelic import frac2, moser_colour, moser_phi, sqrt_2adic
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
WINDOW = {Fr(1, 4), Fr(1, 2), Fr(3, 4)}


def test_the_local_torus_mod_4_is_the_sixth_roots_of_unity():
    """Every norm-one unit a + b omega of Z_2[omega] (a^2 + ab + b^2 = 1) is, mod 4, one of the six
    sixth roots of unity, and (a + 2b)/4 is then 1/4, 1/2 or 3/4.  Mod 2^m for m <= 6 every class
    reduces to one of these six, so nothing deeper changes the value."""
    six = {(1, 0), (0, 1), (3, 1), (3, 0), (0, 3), (1, 3)}     # 1, w, w^2 = w - 1, -1, -w, -w^2 mod 4
    for m in range(2, 7):
        q = 1 << m
        T = [(a, b) for a in range(q) for b in range(q) if (a * a + a * b + b * b - 1) % q == 0]
        assert {(a % 4, b % 4) for a, b in T} == six
    assert {Fr((a + 2 * b) % 4, 4) for a, b in six} == WINDOW


def test_frac2_and_the_square_root_of_33():
    s = sqrt_2adic(33)
    assert (s * s - 33) % (1 << 200) == 0 and s % 4 == 1
    assert frac2(Fr(3, 4)) == Fr(3, 4) and frac2(Fr(5, 3)) == 0 and frac2(Fr(1, 12)) == Fr(3, 4)


def _rotation_units(F):
    r3, r11 = F.sqrt(3), F.sqrt(11)
    q = lambda a, b: F.rational(Fr(a, b))
    R = {"omega": Rotation(q(1, 2), r3 * q(1, 2)), "sigma": Rotation(q(5, 6), r11 * q(1, 6)),
         "lambda": Rotation(q(49, 50), r11 * q(3, 50)), "rho7": Rotation(q(1, 7), r3 * q(4, 7))}
    inv = {k: Rotation(r.cos, -r.sin) for k, r in R.items()}
    u0 = Point(F.one(), F.rational(0))
    units = set()
    for b, c, d in itertools.product(range(-2, 3), range(-1, 2), range(-1, 2)):
        p = u0
        for name, e in (("sigma", b), ("lambda", c), ("rho7", d)):
            for _ in range(abs(e)):
                p = (R[name] if e > 0 else inv[name])(p)
        for _ in range(6):
            units.add(p)
            p = R["omega"](p)
    return units


def test_every_rotation_word_unit_is_in_the_window_at_both_places():
    """Moser's sigma (denominator 6), Exoo-Ismailescu's lambda (25) and rho7 (7): all unit vectors
    they generate have phi in {1/4, 1/2, 3/4}, at both places above 2."""
    F = Field((3, 11))
    units = _rotation_units(F)
    assert len(units) == 270
    for place in (1, -1):
        assert {moser_phi(u, place) for u in units} == WINDOW


def test_the_spindle_and_exoo_ismailescu_H_are_properly_coloured():
    F = Field((3, 11))
    r3, r11 = F.sqrt(3), F.sqrt(11)
    q = lambda a, b: F.rational(Fr(a, b))
    om = Rotation(q(1, 2), r3 * q(1, 2))
    sig = Rotation(q(5, 6), r11 * q(1, 6))
    o = Point(F.rational(0), F.rational(0))
    e = Point(F.one(), F.rational(0))
    a, b = e, om(e)
    tip = a + b
    spindle = [o, a, b, tip] + [sig(p) for p in (a, b, tip)]
    g = build_graph(spindle)
    assert g.m == 11
    col = [moser_colour(p) for p in spindle]
    assert all(col[i] != col[j] for i, j in g.edges())
    d = json.load(open(f"{ROOT}/data/ei_H214.json"))
    Fe = Field(tuple(d["field_generators"]))
    P = [Point(Fe.element([Fr(x, y) for x, y in xy[0]]), Fe.element([Fr(x, y) for x, y in xy[1]]))
         for xy in d["points"]]
    gH = build_graph(P)
    assert gH.m == 1004
    for place in (1, -1):
        colH = [moser_colour(p, place) for p in P]
        assert all(colH[i] != colH[j] for i, j in gH.edges())
