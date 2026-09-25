"""The plane over Q(sqrt3, sqrt11) is 4-colourable, so its chromatic number is exactly 4 (hn/adelic.py).

Moorhouse (2010) left chi(Q(sqrt3, sqrt11)^2) open, Madore (arXiv 1509.07023) proved 4 <= chi <= 5, and
Exoo-Ismailescu (arXiv 1805.00157) asked whether a 5-chromatic unit-distance graph embeds in this plane.
It does not. Write L = Q(sqrt3, sqrt11) and K = L(i). The two places of L over 2 have completion
Q_2(sqrt3) and are inert in K. So every unit vector is a unit of O_w = Z_2[sqrt3][w], with a nonzero
residue in F_4, and the residue of z - rep(z) is a proper 4-colouring. The Moser spindle gives the
lower bound.
"""
import itertools
import json
import os
import random
from fractions import Fraction as Fr

from sympy import QQ, Poly, minimal_polynomial, sqrt, symbols
from sympy.polys.numberfields.primes import prime_decomp

from hn.adelic import PRECISION, q311_colour, q311_coordinates
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
K_SCALE = 64


def _is_square_mod_8(n):
    return n % 8 == 1


def test_the_places_over_2_are_inert_with_residue_field_f4():
    # sqrt33 is 2-adic, so the completion of L over 2 is Q_2(sqrt3), twice.
    assert _is_square_mod_8(33)
    # i lies in Q_2(sqrt3) iff -1 or -3 is a square in Q_2. Neither is (odd squares are 1 mod 8).
    assert not _is_square_mod_8(-1) and not _is_square_mod_8(-3)
    # So K_w = Q_2(sqrt3)(sqrt-3), and -3 = 5 mod 8 makes Q_2(sqrt-3) the unramified quadratic
    # extension. Its residue field is F_2[w] with w^2 + w + 1 = 0, irreducible over F_2.
    assert all((x * x + x + 1) % 2 for x in (0, 1))
    # sqrt3 is a unit (norm -3), and sqrt3 - 1 has norm -2: a uniformiser, so sqrt3 = 1 mod it.
    assert (1 - 3) == -2 and (-3) % 2 == 1


def test_sympy_finds_two_places_over_2_with_residue_field_f2():
    # An independent check of the local facts: in L = Q(sqrt3, sqrt11), 2 = P1^2 P2^2 with residue fields F_2.
    x = symbols("x")
    T = Poly(minimal_polynomial(sqrt(3) + sqrt(11), x), x, domain=QQ)
    assert sorted((P.e, P.f) for P in prime_decomp(2, T)) == [(2, 1), (2, 1)]


def _field_and_units():
    F = Field((3, 11))
    q = lambda a, b: F.rational(Fr(a, b))
    r3, r11, r33 = F.sqrt(3), F.sqrt(11), F.sqrt(33)
    gens = {"w6": Rotation(q(1, 2), r3 * q(1, 2)),        # 60 degrees
            "moser": Rotation(q(5, 6), r11 * q(1, 6)),     # the Moser spindle's rotation
            "345": Rotation(q(3, 5), q(4, 5)),             # not in the Moser field Q(sqrt-3, sqrt-11)
            "s33": Rotation(r33 * q(1, 7), q(4, 7)),       # (sqrt33 + 4i)/7
            "ei": Rotation(q(49, 50), r11 * q(3, 50))}     # Exoo-Ismailescu's lambda
    for r in gens.values():
        assert r.cos * r.cos + r.sin * r.sin == F.one()
    inv = {k: Rotation(r.cos, -r.sin) for k, r in gens.items()}
    e = Point(F.one(), F.rational(0))
    units = set()
    for a, b, c, d in itertools.product(range(-2, 3), range(-1, 2), range(-1, 2), range(-1, 2)):
        p = e
        for name, k in (("moser", a), ("345", b), ("s33", c), ("ei", d)):
            for _ in range(abs(k)):
                p = (gens[name] if k > 0 else inv[name])(p)
        for _ in range(6):
            units.add(p)
            p = gens["w6"](p)
    return F, units


def test_every_unit_vector_is_a_2_adic_unit_with_nonzero_residue():
    F, units = _field_and_units()
    assert len(units) == 810
    for place in (1, -1):
        for u in units:
            co = q311_coordinates(u, place, K_SCALE)
            assert all(v % (1 << K_SCALE) == 0 for v in co), "a unit vector must be 2-adically integral"
            a, b, c, d = ((v >> K_SCALE) & 1 for v in co)
            assert ((a + b) & 1, (c + d) & 1) != (0, 0), "its residue in F_4 must be nonzero"


def test_random_unit_vectors_are_2_adic_units_with_nonzero_residue():
    # Hilbert 90: every unit vector of K = L(i) is t / tbar, i.e. ((a^2 - b^2), 2ab) / (a^2 + b^2).
    F = Field((3, 11))
    basis = (F.one(), F.sqrt(3), F.sqrt(11), F.sqrt(33))
    rng = random.Random(90)
    tested = 0
    while tested < 300:
        a = sum((F.rational(Fr(rng.randint(-9, 9), rng.randint(1, 12))) * e for e in basis), F.zero())
        b = sum((F.rational(Fr(rng.randint(-9, 9), rng.randint(1, 12))) * e for e in basis), F.zero())
        d = a * a + b * b
        if d == F.zero():
            continue
        inv = F.one() / d
        u = Point((a * a - b * b) * inv, F.rational(2) * a * b * inv)
        assert u.x * u.x + u.y * u.y == F.one()
        for place in (1, -1):
            co = q311_coordinates(u, place, K_SCALE)
            assert all(v % (1 << K_SCALE) == 0 for v in co), "a unit vector must be 2-adically integral"
            a0, b0, c0, d0 = ((v >> K_SCALE) & 1 for v in co)
            assert ((a0 + b0) & 1, (c0 + d0) & 1) != (0, 0), "its residue in F_4 must be nonzero"
        tested += 1


def test_the_colouring_is_proper_on_random_unit_steps():
    F, units = _field_and_units()
    units = sorted(units, key=lambda p: (p.fx, p.fy))
    rng = random.Random(311)
    q = lambda a, b: F.rational(Fr(a, b))
    r3, r11, r33 = F.sqrt(3), F.sqrt(11), F.sqrt(33)
    for _ in range(300):
        z = Point(q(rng.randrange(-40, 40), rng.choice([1, 2, 3, 5, 8, 12])) + r33 * q(rng.randrange(-9, 9), 4),
                  r3 * q(rng.randrange(-9, 9), 6) + r11 * q(rng.randrange(-5, 5), rng.choice([1, 3, 16])))
        u = rng.choice(units)
        for place in (1, -1):
            assert q311_colour(z, place) != q311_colour(z + u, place)


def test_the_moser_spindle_needs_four_and_gets_four():
    F = Field((3, 11))
    q = lambda a, b: F.rational(Fr(a, b))
    om = Rotation(q(1, 2), F.sqrt(3) * q(1, 2))
    sig = Rotation(q(5, 6), F.sqrt(11) * q(1, 6))
    o, e = Point(F.rational(0), F.rational(0)), Point(F.one(), F.rational(0))
    a, b = e, om(e)
    spindle = [o, a, b, a + b] + [sig(p) for p in (a, b, a + b)]
    g = build_graph(spindle)
    E = list(g.edges())
    assert len(E) == 11
    assert not any(all(c[i] != c[j] for i, j in E) for c in itertools.product(range(3), repeat=7))
    for place in (1, -1):
        col = [q311_colour(p, place) for p in spindle]
        assert all(col[i] != col[j] for i, j in E)


def test_exoo_ismailescu_H_is_properly_coloured():
    d = json.load(open(os.path.join(ROOT, "data", "ei_H214.json")))
    Fe = Field(tuple(d["field_generators"]))
    P = [Point(Fe.element([Fr(x, y) for x, y in xy[0]]), Fe.element([Fr(x, y) for x, y in xy[1]]))
         for xy in d["points"]]
    g = build_graph(P)
    assert g.m == 1004
    for place in (1, -1):
        col = [q311_colour(p, place) for p in P]
        assert all(col[i] != col[j] for i, j in g.edges())
        assert len(set(col)) == 4


def test_pairs_at_8_over_3_are_alike_and_at_sqrt_11_over_3_apart():
    """Exoo-Ismailescu's G_40 forces a pair at distance 8/3 alike in every 4-colouring with no
    monochromatic pair at distance sqrt(11/3), and Parts (Polymath16, July 2019) chained such pairs into
    alike pairs at every distance 8/9^n, which sum to 1. His limit argument is not a proof. The colouring
    here has no monochromatic pair at sqrt(11/3) = sqrt33/3, a 2-adic unit, so E-I's hypothesis holds;
    it makes EVERY pair at distance 8/9^n alike, since (8/9^n) u lies in 8 O_w, and it is still proper."""
    F, units = _field_and_units()
    units = sorted(units, key=lambda p: (p.fx, p.fy))
    rng = random.Random(83)
    q = lambda a, b: F.rational(Fr(a, b))
    for _ in range(200):
        z = Point(q(rng.randrange(-40, 40), rng.choice([1, 2, 3, 7])) + F.sqrt(11) * q(rng.randrange(-5, 5), 6),
                  F.sqrt(3) * q(rng.randrange(-9, 9), 4) + F.sqrt(33) * q(rng.randrange(-3, 3), 9))
        u = rng.choice(units)
        n = rng.randrange(0, 4)
        for place in (1, -1):
            c = q311_colour(z, place)
            assert q311_colour(z + u.scaled(q(8, 9 ** n)), place) == c
            assert q311_colour(z + u.scaled(F.sqrt(33) * q(1, 3)), place) != c
    assert PRECISION > 2 * K_SCALE
