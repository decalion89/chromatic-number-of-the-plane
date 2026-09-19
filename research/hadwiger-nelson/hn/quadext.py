"""K(sqrt(d)) over a cyclotomic field: the smallest home for an order-5 spindle.

Two requirements pull in opposite directions, and only one field satisfies
both cheaply.

*The radius.*  Targets have to sit on the order-n magic circle, whose squared
radius 1/(2 - zeta_n - zeta_n^-1) lives in Q(zeta_n)^+, of degree phi(n)/2.
The points reaching it therefore need a field of at least that degree, and
degree is what makes a group high-rank -- which makes its balls all boundary,
with nothing left after peeling.  Only orders 3 and 5 have their radius in a
degree-2 real field: 1 and 2 respectively, against 3 at order 7, 5 at order 11
and 6 at order 13.  Order 3 is the classical wall, capacity 2.  **Order 5 is
the only reachable escape, and it has capacity 3.**

*The structure.*  A Moser spindle needs omega for its rhombi and the rotation
(5 + sqrt(-11))/6 to close them, and Q(zeta_15) -- which has zeta_5, sqrt5 and
omega -- has no sqrt(-11).  Its walks come out 3-chromatic however many unit
steps are thrown at them, 60000 vertices included.

But only sqrt(-11) is missing, not all of Q(zeta_11).  Adjoining it to
Q(zeta_15) gives degree 16, not the 80 of Q(zeta_165), and carries everything:
omega for the triangles, zeta_5 for the order-5 rotations, sqrt5 for the
order-5 radius, sqrt(-11) for the spindle.

Elements are flat tuples, the base field's coefficients followed by the
coefficients of the sqrt(d) part, so they work unchanged with the integer walk
and the graph builders.
"""
from __future__ import annotations

from fractions import Fraction
from math import sqrt
from typing import List, Tuple

from .cyclotomic import CycloField


class QuadExtField:
    """Q(zeta_n)(sqrt(d)) with d a negative integer, so sqrt(d) is imaginary."""

    def __init__(self, base: CycloField, d: int):
        if d >= 0:
            raise ValueError("d must be negative so sqrt(d) is imaginary")
        self.base = base
        self.d = d
        self.half = base.degree
        self.degree = 2 * base.degree
        self._sqrt_d = complex(0.0, sqrt(-d))

    # -- splitting and joining ----------------------------------------------
    def _split(self, x) -> Tuple[Tuple, Tuple]:
        return x[: self.half], x[self.half:]

    def _join(self, a, b) -> Tuple:
        return tuple(a) + tuple(b)

    # -- arithmetic ----------------------------------------------------------
    def zero(self) -> Tuple:
        return (Fraction(0),) * self.degree

    def one(self) -> Tuple:
        return self._join(self.base.one(), self.base.zero())

    def rational(self, q) -> Tuple:
        return self._join(self.base.rational(q), self.base.zero())

    def embed(self, a) -> Tuple:
        """A base-field element, as an element here."""
        return self._join(a, self.base.zero())

    def radical(self) -> Tuple:
        """sqrt(d) itself."""
        return self._join(self.base.zero(), self.base.one())

    def add(self, x, y) -> Tuple:
        return tuple(a + b for a, b in zip(x, y))

    def sub(self, x, y) -> Tuple:
        return tuple(a - b for a, b in zip(x, y))

    def neg(self, x) -> Tuple:
        return tuple(-a for a in x)

    def scale(self, x, q) -> Tuple:
        q = Fraction(q)
        return tuple(a * q for a in x)

    def mul(self, x, y) -> Tuple:
        a, b = self._split(x)
        c, e = self._split(y)
        B = self.base
        # (a + b r)(c + e r) = (a c + d b e) + (a e + b c) r,  r^2 = d
        real = B.add(B.mul(a, c), B.scale(B.mul(b, e), self.d))
        rad = B.add(B.mul(a, e), B.mul(b, c))
        return self._join(real, rad)

    def conj(self, x) -> Tuple:
        """Complex conjugation: it flips sqrt(d), which is purely imaginary."""
        a, b = self._split(x)
        B = self.base
        return self._join(B.conj(a), B.neg(B.conj(b)))

    def norm2(self, x) -> Tuple:
        return self.mul(x, self.conj(x))

    def is_unit_apart(self, x, y) -> bool:
        return self.norm2(self.sub(x, y)) == self.one()

    def to_complex(self, x) -> complex:
        a, b = self._split(x)
        return self.base.to_complex(a) + self.base.to_complex(b) * self._sqrt_d

    # -- the elements a construction asks for --------------------------------
    def zeta(self, k: int = 1) -> Tuple:
        return self.embed(self.base.zeta(k))

    def root_of_unity(self, order: int, k: int = 1) -> Tuple:
        if self.base.n % order:
            raise ValueError(f"zeta_{order} is not in Q(zeta_{self.base.n})")
        return self.embed(self.base.zeta((self.base.n // order) * k))

    def moser_rotation(self) -> Tuple:
        """(5 + sqrt(-11))/6, which needs d = -11."""
        if self.d != -11:
            raise ValueError("the Moser rotation needs sqrt(-11)")
        return self.scale(self.add(self.rational(5), self.radical()),
                          Fraction(1, 6))

    def unit_steps(self, rot=None, m_max: int = 0) -> List[Tuple]:
        """Roots of unity of the base, times powers of a rotation."""
        base_steps = [self.embed(tuple(Fraction(c) for c in u))
                      for u in self.base.unit_steps()]
        if rot is None or m_max == 0:
            return base_steps
        powers = [self.one()]
        for _ in range(m_max):
            powers.append(self.mul(powers[-1], rot))
        inv = self.conj(rot)
        for _ in range(m_max):
            powers.append(self.mul(powers[-1], inv))
        seen = {}
        for m in powers:
            for s in base_steps:
                z = self.mul(s, m)
                seen[z] = None
                seen[self.neg(z)] = None
        return [z for z in seen if self.norm2(z) == self.one()]


def degrey_field() -> "QuadExtField":
    """Q(zeta_15, sqrt(-7), sqrt(-11)): de Grey's construction and zeta_15 at once.

    Read as pairs (cos, sin) his rotations look like they need the real field
    Q(sqrt3, sqrt5, sqrt7, sqrt11).  Read as complex numbers they are

        2 arcsin(1/4)       -> (7 + sqrt(-15)) / 8
        2 arcsin(1/8)       -> (31 + 3 sqrt(-7)) / 32
        pi/2 - arcsin(1/8)  -> (-1 + 3 sqrt(-7)) / 8

    each of modulus exactly 1, and the Moser rotation (5 + sqrt(-11))/6 with
    them.  The radicals are sqrt(-15), sqrt(-7) and sqrt(-11), and sqrt(-15) =
    sqrt(-3) sqrt(5) is already inside Q(zeta_15).

    So two quadratic steps over the cyclotomic base carry the whole 2018
    construction -- in a field that also has zeta_15, whose magic circle has
    capacity 4 against the multiquadratic 2.  The multiquadratic field is not
    where this lives; it is only where the coordinates were written down.
    """
    return QuadExtField(QuadExtField(CycloField(15), -7), -11)


def degrey_rotations(L: "QuadExtField") -> dict:
    """The four rotations, as elements of `degrey_field()`."""
    inner = L.base
    base = inner.base
    r7 = L.embed(inner.radical())
    r11 = L.radical()
    r3 = L.embed(inner.embed(base.sqrt_disc(3)))      # sqrt(-3), 3 = 3 mod 4
    r5 = L.embed(inner.embed(base.sqrt_disc(5)))      # sqrt(5),  5 = 1 mod 4
    r15 = L.mul(r3, r5)                               # sqrt(-15)
    return {
        "spindle_2": L.scale(L.add(L.rational(7), r15), Fraction(1, 8)),
        "spindle_4": L.scale(L.add(L.rational(31), L.scale(r7, 3)), Fraction(1, 32)),
        "quarter": L.scale(L.add(L.rational(-1), L.scale(r7, 3)), Fraction(1, 8)),
        "moser": L.scale(L.add(L.rational(5), r11), Fraction(1, 6)),
    }
