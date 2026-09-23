"""The three- and six-fold disjunctive spindles, checked in exact arithmetic.

De Grey's spindle consumes a pair FORCED equal in every k-colouring.  These two
consume something strictly weaker -- a disjunction

        c(v) = c(q1)   OR   c(v) = c(q2)

-- at the price of more copies.  Everything below is the geometry and the
combinatorics of that trade, recomputed rather than quoted: the distances are
exact field elements compared to 1, and the pigeonhole is a brute-force sweep
over every way the copies could choose their options.

The argument.  Put s rotations about v, by angles theta_0 .. theta_{s-1}, and
let U be the union of the s rotated copies of H.  Every copy contains v, which
every rotation fixes, so copy j gives c(v) = c(rho_j(q_i)) for some option
i(j).  Two copies holding the SAME option collide when their two points are a
unit apart, and since both lie on the circle of radius r_i about v,

        |rho_j(q_i) - rho_k(q_i)| = 2 r_i |sin((theta_j - theta_k)/2)|

which is 1 exactly when theta_j - theta_k = +- alpha_i, the spindle angle
2 arcsin(1/(2 r_i)) of that radius.  So option i's copies must form an
independent set in the graph G_i on {0..s-1} joining angles that differ by
+- alpha_i, and the argument closes when the options cannot cover all s copies:

        sum_i  alpha(G_i)  <  s

Each G_i has maximum degree 2, so alpha(G_i) >= s/3 with equality only for
disjoint triangles, and >= s/2 for a union of paths.  Two unions of paths give
at least s and never close, so at least one G_i must contain an odd cycle --
its alpha must be a rational multiple of 2 pi.  The triangle needs
alpha = 2 pi / 3, hence r = 1/sqrt3, which every field in this project holds;
the next odd cycle is the pentagon, needing sin 72 = sqrt(10 + 2 sqrt5)/4, a
nested radical no multiquadratic field contains.  1/sqrt3 is the only special
radius available here, and at least one endpoint must sit on it.

Three copies:  angles {0, 2pi/3, 4pi/3} and BOTH endpoints at 1/sqrt3.
               G_1 = G_2 = the triangle, alpha 1 each, 1 + 1 < 3.

Six copies:    angles {0, 2pi/3, 4pi/3} u {beta, beta+2pi/3, beta+4pi/3},
               q1 at 1/sqrt3 and q2 at ANY radius > 1/2.
               G_1 = two triangles (alpha 2), G_2 = a perfect matching
               (alpha 3), and 2 + 3 < 6.

The six-copy version is the one worth having: it asks for one special distance
instead of two.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from itertools import product

import pytest

from hn.field import Field
from hn.geometry import Point, Rotation, rotation_joining

FIELD = Field((3, 5))
ONE = FIELD.rational(Fr(1))
ZERO = FIELD.rational(Fr(0))
ORIGIN = Point(ZERO, ZERO)

# |v - q1|^2 = 1/3 exactly:  q1 = (sqrt3/3, 0)
Q1 = Point(FIELD.sqrt(3) * FIELD.rational(Fr(1, 3)), ZERO)
# |v - q2|^2 = 4, so beta = 2 arcsin(1/4), cos 7/8, sin sqrt15/8 -- not a
# multiple of 60 degrees, so the six rotations really are six.
Q2 = Point(FIELD.rational(Fr(2)), ZERO)

ROT120 = Rotation(FIELD.rational(Fr(-1, 2)), FIELD.sqrt(3) * FIELD.rational(Fr(1, 2)))
ROTBETA = rotation_joining(Fr(4), FIELD)


def _compose(a: Rotation, b: Rotation) -> Rotation:
    return Rotation(a.cos * b.cos - a.sin * b.sin, a.cos * b.sin + a.sin * b.cos)


def _identity() -> Rotation:
    return Rotation(ONE, ZERO)


def _cosets():
    """The six rotations, as (coset, power) -> rotation about the origin."""
    out, cur = {}, _identity()
    triangle = []
    for _ in range(3):
        triangle.append(cur)
        cur = _compose(cur, ROT120)
    for c, shift in enumerate((_identity(), ROTBETA)):
        for p, t in enumerate(triangle):
            out[(c, p)] = _compose(shift, t)
    return out


def test_q1_sits_at_one_third_and_q2_above_a_quarter():
    assert (Q1 - ORIGIN).norm2() == FIELD.rational(Fr(1, 3))
    # the spindle angle of r2 exists only for r2 > 1/2
    assert float((Q2 - ORIGIN).norm2()) > 0.25


def test_the_six_rotations_are_distinct():
    rots = _cosets()
    seen = {(r.cos, r.sin) for r in rots.values()}
    assert len(seen) == 6


def test_option_one_collides_exactly_within_a_coset():
    """q1 at 1/sqrt3: its clash graph is the two triangles, nothing else.

    Two copies clash on q1 when their images of q1 are a unit apart.  The
    triangle inscribed in the circle of radius r has side r*sqrt3, so at
    r = 1/sqrt3 that side is exactly 1 -- and the coset structure means that
    happens for the two copies of the SAME coset and no others.
    """
    rots = _cosets()
    keys = sorted(rots)
    edges = set()
    for a in keys:
        for b in keys:
            if a >= b:
                continue
            d = (rots[a](Q1) - rots[b](Q1)).norm2()
            if d == ONE:
                edges.add((a, b))
    expected = {((c, p), (c, q)) for c in (0, 1)
                for p in range(3) for q in range(3) if p < q}
    assert edges == expected
    assert len(edges) == 6          # two vertex-disjoint triangles


def test_option_two_collides_exactly_across_the_cosets():
    """q2 at radius 2: its clash graph is the perfect matching.

    beta is the spindle angle of r2, so the two copies whose angles differ by
    exactly beta -- one from each coset, at the same power of the 120 degree
    rotation -- put their images of q2 a unit apart.
    """
    rots = _cosets()
    keys = sorted(rots)
    edges = set()
    for a in keys:
        for b in keys:
            if a >= b:
                continue
            if (rots[a](Q2) - rots[b](Q2)).norm2() == ONE:
                edges.add((a, b))
    expected = {((0, p), (1, p)) for p in range(3)}
    assert edges == expected
    assert len(edges) == 3          # a perfect matching on the six copies


def _independent(chosen, edges):
    s = set(chosen)
    return not any(a in s and b in s for a, b in edges)


def test_two_options_cannot_cover_six_copies():
    """The pigeonhole, by brute force over all 64 ways to choose.

    Every copy must take one of the two options.  Option 1's copies have to be
    independent in the two triangles (at most 2 of them) and option 2's in the
    matching (at most 3), so at most 5 of the 6 copies can be covered.  The
    sweep confirms there is no exception.
    """
    rots = _cosets()
    keys = sorted(rots)
    g1 = {((c, p), (c, q)) for c in (0, 1)
          for p in range(3) for q in range(3) if p < q}
    g2 = {((0, p), (1, p)) for p in range(3)}
    survivors = 0
    for choice in product((1, 2), repeat=6):
        a1 = [k for k, ch in zip(keys, choice) if ch == 1]
        a2 = [k for k, ch in zip(keys, choice) if ch == 2]
        if _independent(a1, g1) and _independent(a2, g2):
            survivors += 1
    assert survivors == 0
    # and the bound is tight: drop one copy and an assignment appears
    loose = 0
    for choice in product((1, 2), repeat=5):
        a1 = [k for k, ch in zip(keys[:5], choice) if ch == 1]
        a2 = [k for k, ch in zip(keys[:5], choice) if ch == 2]
        if _independent(a1, g1) and _independent(a2, g2):
            loose += 1
    assert loose > 0


def test_three_copies_suffice_when_both_endpoints_are_special():
    """The cheaper version: both at 1/sqrt3, three copies, one triangle each."""
    q2 = ROT120(Q1)                          # another point at 1/sqrt3
    assert (q2 - ORIGIN).norm2() == FIELD.rational(Fr(1, 3))
    rots, cur = [], _identity()
    for _ in range(3):
        rots.append(cur)
        cur = _compose(cur, ROT120)
    for q in (Q1, q2):
        for a in range(3):
            for b in range(a + 1, 3):
                assert (rots[a](q) - rots[b](q)).norm2() == ONE
    # every pair of the three copies clashes on either option, so two options
    # can cover at most 1 + 1 = 2 < 3
    survivors = [c for c in product((1, 2), repeat=3)
                 if c.count(1) <= 1 and c.count(2) <= 1]
    assert survivors == []


@pytest.mark.parametrize("r2_squared", [1, 3, 4, 16])
def test_the_second_radius_is_genuinely_free(r2_squared):
    """Any r2 whose spindle angle the field holds works -- that is the point.

    The matching is built from beta = 2 arcsin(1/(2 r2)) whatever r2 is, so the
    second endpoint carries no distance constraint beyond r2 > 1/2 and the one
    radical sqrt(4 r2^2 - 1) that every spindle letter needs.  The four radii
    here are the letters of de Grey's own alphabet: 3, 15 = 3*5, 63 = 9*7 and
    11.  (At r2^2 = 1 the two cosets merge into C6 and the matching becomes the
    long diagonals -- still a perfect matching, still independence 3.)
    """
    field = Field((3, 5, 7, 11))
    zero = field.zero()
    one = field.rational(Fr(1))
    rot120 = Rotation(field.rational(Fr(-1, 2)),
                      field.sqrt(3) * field.rational(Fr(1, 2)))
    beta = rotation_joining(Fr(r2_squared), field)
    q2 = Point(field.sqrt(r2_squared), zero)
    assert (q2 - Point(zero, zero)).norm2() == field.rational(Fr(r2_squared))
    tri, cur = [], Rotation(one, zero)
    for _ in range(3):
        tri.append(cur)
        cur = Rotation(cur.cos * rot120.cos - cur.sin * rot120.sin,
                       cur.cos * rot120.sin + cur.sin * rot120.cos)
    for t in tri:
        b = Rotation(beta.cos * t.cos - beta.sin * t.sin,
                     beta.cos * t.sin + beta.sin * t.cos)
        assert (t(q2) - b(q2)).norm2() == one
