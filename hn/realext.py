"""A real quadratic extension of a multiquadratic field.

`hn.field.Field` is multiquadratic by construction, and multiquadratic real
fields are totally real: every Galois conjugate of every element is real. That
is a genuine restriction, and it is what stopped the four-point configuration.

Its pivot sits at squared distance (7 + sqrt(33))/6 from the second leg, and
the rotation closing that circle has

    sin^2 = (-66 + 30 sqrt(33)) / 256,

whose conjugate, with sqrt(33) negated, is negative. A square root of a
non-totally-positive element lies in no totally real field, so no amount of
widening inside `Field` will ever name it -- while the rotation itself plainly
exists in the plane, at about 40.1 degrees.

Adjoining one real square root fixes it. Elements are pairs (a, b) standing
for a + b sqrt(v) with v a positive element of the base, arithmetic is the
obvious thing, and inversion is the conjugate trick
1/(a + b r) = (a - b r)/(a^2 - v b^2). Nothing here is complex: the extension
is real, so the points it builds are ordinary plane points and every distance
comparison stays exact.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Optional

from .field import Field, FieldElement


class RealQuadExt:
    """K(sqrt(v)) for K multiquadratic and v a positive element of K."""

    __slots__ = ("base", "v", "_vf")

    def __init__(self, base: Field, v: FieldElement):
        if v.field is not base:
            raise TypeError("v must live in the base field")
        self.base = base
        self.v = v
        self._vf = float(v)
        if self._vf <= 0:
            raise ValueError("v must be positive in the real embedding")

    def __eq__(self, other) -> bool:
        return (isinstance(other, RealQuadExt) and other.base is self.base
                and other.v == self.v)

    def __hash__(self) -> int:
        return hash(("RealQuadExt", self.v))

    def __repr__(self) -> str:
        return f"{self.base}(sqrt({self.v}))"

    # -- constructors --------------------------------------------------------
    def zero(self) -> "RealExtElement":
        return RealExtElement(self, self.base.rational(0), self.base.rational(0))

    def rational(self, q) -> "RealExtElement":
        return RealExtElement(self, self.base.rational(q), self.base.rational(0))

    def embed(self, a: FieldElement) -> "RealExtElement":
        return RealExtElement(self, a, self.base.rational(0))

    def radical(self) -> "RealExtElement":
        """sqrt(v) itself."""
        return RealExtElement(self, self.base.rational(0), self.base.rational(1))

    def sqrt(self, d: int) -> "RealExtElement":
        return self.embed(self.base.sqrt(d))


class RealExtElement:
    """a + b*sqrt(v), with a and b in the base field."""

    __slots__ = ("field", "a", "b")

    def __init__(self, field: RealQuadExt, a: FieldElement, b: FieldElement):
        self.field = field
        self.a = a
        self.b = b

    # -- helpers -------------------------------------------------------------
    def _coerce(self, other):
        if isinstance(other, RealExtElement):
            if other.field != self.field:
                raise TypeError("mixing different extensions")
            return other
        if isinstance(other, (int, Fraction)):
            return self.field.rational(other)
        if isinstance(other, FieldElement):
            return self.field.embed(other)
        return NotImplemented

    def is_rational(self) -> bool:
        return self.b == 0 and self.a.is_rational()

    @property
    def c(self):
        """The base coefficients, so rational elements read like Field ones."""
        return self.a.c

    # -- arithmetic ----------------------------------------------------------
    def __add__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return RealExtElement(self.field, self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return RealExtElement(self.field, -self.a, -self.b)

    def __sub__(self, other):
        o = self._coerce(other)
        return NotImplemented if o is NotImplemented else self + (-o)

    def __rsub__(self, other):
        o = self._coerce(other)
        return NotImplemented if o is NotImplemented else o + (-self)

    def __mul__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return RealExtElement(self.field,
                              self.a * o.a + self.field.v * self.b * o.b,
                              self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def inverse(self) -> "RealExtElement":
        """(a - b r) / (a^2 - v b^2); the norm is zero only for zero itself."""
        norm = self.a * self.a - self.field.v * self.b * self.b
        if norm == 0:
            raise ZeroDivisionError("element has norm zero")
        inv = norm.inverse()
        return RealExtElement(self.field, self.a * inv, -self.b * inv)

    def __truediv__(self, other):
        o = self._coerce(other)
        return NotImplemented if o is NotImplemented else self * o.inverse()

    def __rtruediv__(self, other):
        o = self._coerce(other)
        return NotImplemented if o is NotImplemented else o * self.inverse()

    # -- comparison ----------------------------------------------------------
    def __eq__(self, other) -> bool:
        if isinstance(other, (int, Fraction)):
            return self.b == 0 and self.a == other
        if isinstance(other, FieldElement):
            return self.b == 0 and self.a == other
        return (isinstance(other, RealExtElement) and self.field == other.field
                and self.a == other.a and self.b == other.b)

    def __hash__(self) -> int:
        return hash((self.a, self.b))

    def __float__(self) -> float:
        return float(self.a) + float(self.b) * self.field._vf ** 0.5

    def __repr__(self) -> str:
        if self.b == 0:
            return repr(self.a)
        return f"({self.a}) + ({self.b})*sqrt({self.field.v})"


def extend_for_sqrt(v: FieldElement) -> Optional[RealQuadExt]:
    """The extension holding sqrt(v), or None when v is not positive."""
    if float(v) <= 0:
        return None
    return RealQuadExt(v.field, v)
