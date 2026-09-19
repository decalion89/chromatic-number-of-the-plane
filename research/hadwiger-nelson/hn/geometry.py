"""Points and rotations of the plane, with exact coordinates.

A `Point` has both coordinates in one `Field`, so `dist2` is exact and the
predicate "these two points are at distance exactly 1" is decidable.

A `Rotation` is a pair (cos, sin) of field elements with cos^2 + sin^2 = 1,
which is checked on construction -- a rotation that is only approximately a
rotation would silently corrupt every edge downstream.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable, Sequence

from .field import Field, FieldElement, QSQRT3_11

__all__ = ["Point", "Rotation", "origin", "eisenstein", "ROT60", "SPINDLE", "rotation_joining"]


class Point:
    """A point of the plane with exact coordinates."""

    __slots__ = ("x", "y", "_hash", "_fx", "_fy")

    def __init__(self, x: FieldElement, y: FieldElement):
        if x.field != y.field:
            raise TypeError("coordinates live in different fields")
        self.x = x
        self.y = y
        self._hash = None
        self._fx = None
        self._fy = None

    @property
    def field(self) -> Field:
        return self.x.field

    @property
    def fx(self) -> float:
        if self._fx is None:
            self._fx = float(self.x)
        return self._fx

    @property
    def fy(self) -> float:
        if self._fy is None:
            self._fy = float(self.y)
        return self._fy

    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y)

    def __neg__(self) -> "Point":
        return Point(-self.x, -self.y)

    def scaled(self, k) -> "Point":
        f = self.field.rational(k) if isinstance(k, (int, Fraction)) else k
        return Point(self.x * f, self.y * f)

    def dist2(self, other: "Point") -> FieldElement:
        dx = self.x - other.x
        dy = self.y - other.y
        return dx * dx + dy * dy

    def is_unit_apart(self, other: "Point") -> bool:
        """Exact test.  This is the only place an edge is ever decided."""
        return self.dist2(other) == 1

    def norm2(self) -> FieldElement:
        return self.x * self.x + self.y * self.y

    def __eq__(self, other) -> bool:
        return isinstance(other, Point) and self.x == other.x and self.y == other.y

    def __hash__(self) -> int:
        if self._hash is None:
            self._hash = hash((self.x, self.y))
        return self._hash

    def __repr__(self) -> str:
        return f"({self.x}, {self.y})"

    def approx(self) -> tuple:
        return (self.fx, self.fy)


class Rotation:
    """Rotation about the origin by an angle with cos, sin in the field."""

    __slots__ = ("cos", "sin")

    def __init__(self, cos: FieldElement, sin: FieldElement, check: bool = True):
        if check and not (cos * cos + sin * sin == 1):
            raise ValueError(f"cos^2+sin^2 = {cos*cos+sin*sin}, not 1 -- not a rotation")
        self.cos = cos
        self.sin = sin

    @property
    def field(self) -> Field:
        return self.cos.field

    def __call__(self, p: Point) -> Point:
        return Point(self.cos * p.x - self.sin * p.y, self.sin * p.x + self.cos * p.y)

    def about(self, pivot: Point):
        """The same rotation, but centred on `pivot`."""

        def rotate(p: Point, _r=self, _c=pivot) -> Point:
            return _r(p - _c) + _c

        return rotate

    def inverse(self) -> "Rotation":
        return Rotation(self.cos, -self.sin, check=False)

    def __mul__(self, other: "Rotation") -> "Rotation":
        """Composition: angles add."""
        return Rotation(
            self.cos * other.cos - self.sin * other.sin,
            self.sin * other.cos + self.cos * other.sin,
            check=False,
        )

    def __pow__(self, n: int) -> "Rotation":
        if n < 0:
            return self.inverse() ** (-n)
        f = self.field
        result = Rotation(f.one(), f.zero(), check=False)
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def __eq__(self, other) -> bool:
        return isinstance(other, Rotation) and self.cos == other.cos and self.sin == other.sin

    def __hash__(self) -> int:
        return hash((self.cos, self.sin))

    def __repr__(self) -> str:
        import math

        return f"Rotation({math.degrees(math.atan2(float(self.sin), float(self.cos))):.4f} deg)"


def origin(field: Field = QSQRT3_11) -> Point:
    return Point(field.zero(), field.zero())


def eisenstein(a: int, b: int, field: Field = QSQRT3_11) -> Point:
    """The triangular-lattice point a + b*omega, omega = exp(i*pi/3).

    Nearest-neighbour distance is 1, so the unit-distance graph on the lattice
    is the triangular lattice itself (chromatic number 3).
    """
    half = field.rational(Fraction(1, 2))
    return Point(field.rational(a) + field.rational(b) * half, field.sqrt(3) * half * field.rational(b))


def _rot60(field: Field = QSQRT3_11) -> Rotation:
    half = field.rational(Fraction(1, 2))
    return Rotation(half, field.sqrt(3) * half)


def _spindle(field: Field = QSQRT3_11) -> Rotation:
    """Rotation by arccos(5/6).

    Two points at distance sqrt(3) from the pivot land at distance exactly 1
    from each other:  3*(2 - 2*(5/6)) = 1.  This is the hinge of the Moser
    spindle and of every 5-chromatic graph built since.
    """
    sixth = field.rational(Fraction(1, 6))
    return Rotation(field.rational(Fraction(5, 6)), field.sqrt(11) * sixth)


ROT60 = _rot60()
SPINDLE = _spindle()


def rotation_joining(d2, field: Field = QSQRT3_11) -> Rotation:
    """The rotation that brings two points at squared distance `d2` from the
    pivot to distance exactly 1 apart.

    From |r*e^{i t} - r|^2 = 2 r^2 (1 - cos t) = 1 we get cos t = 1 - 1/(2 d2)
    and sin t = sqrt(1 - cos^2 t), which must be representable in the field.
    """
    d2 = Fraction(d2)
    c = Fraction(1) - Fraction(1, 2) / d2
    s2 = 1 - c * c
    num, den = s2.numerator, s2.denominator
    # sin = sqrt(num/den) = sqrt(num*den)/den
    root = field.sqrt(num * den)
    return Rotation(field.rational(c), root * field.rational(Fraction(1, den)))
