"""Exact arithmetic in multiquadratic number fields Q(sqrt(d_1), ..., sqrt(d_k)).

Every coordinate of every vertex we ever construct lives in such a field, so
distances can be compared to 1 *exactly*.  Floating point is used only as a
pre-filter; it never decides whether an edge exists.

An element is a vector of 2**k rationals over the basis

    { sqrt(prod_{i in S} d_i) : S subset of {0..k-1} }

indexed by the bitmask of S.  The basis is a genuine Q-basis as long as the
generators are squarefree and pairwise coprime, which `Field` enforces.

Multiplication is driven by  sqrt(A) * sqrt(B) = prod_{i in A&B} d_i * sqrt(A^B),
i.e. the exponent sets combine by symmetric difference.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt
from typing import Iterable, Sequence

__all__ = ["Field", "FieldElement", "QSQRT3_11", "embed"]


def _is_squarefree(n: int) -> bool:
    if n < 1:
        return False
    d = 2
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 1
    return True


class Field:
    """A multiquadratic field, and the factory for its elements."""

    __slots__ = ("gens", "rank", "dim", "_prod", "_names")

    def __init__(self, gens: Sequence[int]):
        gens = tuple(int(g) for g in gens)
        for g in gens:
            if g < 2 or not _is_squarefree(g):
                raise ValueError(f"generator {g} must be squarefree and >= 2")
        for i, a in enumerate(gens):
            for b in gens[i + 1 :]:
                from math import gcd

                if gcd(a, b) != 1:
                    raise ValueError(f"generators {a} and {b} must be coprime")
        self.gens = gens
        self.rank = len(gens)
        self.dim = 1 << self.rank
        # _prod[m] = product of the generators selected by bitmask m
        prod = [1] * self.dim
        for m in range(self.dim):
            p = 1
            for i in range(self.rank):
                if m >> i & 1:
                    p *= gens[i]
            prod[m] = p
        self._prod = tuple(prod)
        names = []
        for m in range(self.dim):
            names.append("1" if m == 0 else f"sqrt({self._prod[m]})")
        self._names = tuple(names)

    # -- constructors ----------------------------------------------------
    def zero(self) -> "FieldElement":
        return FieldElement(self, (Fraction(0),) * self.dim)

    def one(self) -> "FieldElement":
        return self.rational(1)

    def rational(self, value) -> "FieldElement":
        c = [Fraction(0)] * self.dim
        c[0] = Fraction(value)
        return FieldElement(self, tuple(c))

    def sqrt(self, d: int) -> "FieldElement":
        """The element sqrt(d) for any d that is a product of distinct generators."""
        for m in range(self.dim):
            if self._prod[m] == d:
                c = [Fraction(0)] * self.dim
                c[m] = Fraction(1)
                return FieldElement(self, tuple(c))
        # allow a square times a representable radical, e.g. sqrt(12) = 2 sqrt(3)
        r = isqrt(d)
        while r > 1:
            if d % (r * r) == 0:
                inner = d // (r * r)
                for m in range(self.dim):
                    if self._prod[m] == inner:
                        c = [Fraction(0)] * self.dim
                        c[m] = Fraction(r)
                        return FieldElement(self, tuple(c))
            r -= 1
        raise ValueError(f"sqrt({d}) is not in {self}")

    def element(self, coeffs: Iterable) -> "FieldElement":
        c = [Fraction(x) for x in coeffs]
        if len(c) != self.dim:
            raise ValueError(f"expected {self.dim} coefficients, got {len(c)}")
        return FieldElement(self, tuple(c))

    def extend(self, *extra: int) -> "Field":
        """The smallest multiquadratic field containing this one and sqrt(extra)."""
        return Field(self.gens + tuple(g for g in extra if g not in self.gens))

    def __repr__(self) -> str:
        return "Q(" + ", ".join(f"sqrt({g})" for g in self.gens) + ")"

    def __eq__(self, other) -> bool:
        return isinstance(other, Field) and self.gens == other.gens

    def __hash__(self) -> int:
        return hash(("Field", self.gens))


class FieldElement:
    """An exact element of a `Field`.  Immutable and hashable."""

    __slots__ = ("field", "c", "_hash")

    def __init__(self, field: Field, coeffs: tuple):
        self.field = field
        self.c = coeffs
        self._hash = None

    # -- helpers ---------------------------------------------------------
    def _coerce(self, other):
        if isinstance(other, FieldElement):
            if other.field != self.field:
                raise TypeError(f"mixing {self.field} and {other.field}")
            return other
        if isinstance(other, (int, Fraction)):
            return self.field.rational(other)
        return NotImplemented

    # -- arithmetic ------------------------------------------------------
    def __add__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return FieldElement(self.field, tuple(a + b for a, b in zip(self.c, o.c)))

    __radd__ = __add__

    def __neg__(self):
        return FieldElement(self.field, tuple(-a for a in self.c))

    def __sub__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return FieldElement(self.field, tuple(a - b for a, b in zip(self.c, o.c)))

    def __rsub__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return o - self

    def __mul__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        f = self.field
        prod = f._prod
        out = [Fraction(0)] * f.dim
        sc, oc = self.c, o.c
        for i, a in enumerate(sc):
            if not a:
                continue
            for j, b in enumerate(oc):
                if not b:
                    continue
                out[i ^ j] += a * b * prod[i & j]
        return FieldElement(f, tuple(out))

    __rmul__ = __mul__

    def inverse(self) -> "FieldElement":
        """Multiplicative inverse.

        The Galois group of Q(sqrt(d_1),...,sqrt(d_k))/Q is (Z/2)^k, acting by
        independent sign flips of the generators.  Multiplying by every
        non-identity conjugate lands in Q (that product is the norm), so
        x^-1 = (prod of conjugates) / Norm(x).
        """
        f = self.field
        acc = f.one()
        for sign_mask in range(1, f.dim):
            acc = acc * self.conjugate(sign_mask)
        norm = (self * acc).c
        if any(norm[m] for m in range(1, f.dim)):
            raise ArithmeticError("norm is not rational; field inconsistency")
        n = norm[0]
        if n == 0:
            raise ZeroDivisionError("element is zero")
        return acc * f.rational(Fraction(1, 1) / n)

    def conjugate(self, gen_mask: int) -> "FieldElement":
        """Flip the sign of sqrt(d_i) for every i selected by `gen_mask`."""
        out = []
        for m, a in enumerate(self.c):
            out.append(-a if bin(m & gen_mask).count("1") % 2 else a)
        return FieldElement(self.field, tuple(out))

    def __truediv__(self, other):
        o = self._coerce(other)
        if o is NotImplemented:
            return NotImplemented
        return self * o.inverse()

    def __pow__(self, n: int):
        if n < 0:
            return self.inverse() ** (-n)
        result = self.field.one()
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    # -- predicates ------------------------------------------------------
    def is_zero(self) -> bool:
        return not any(self.c)

    def is_rational(self) -> bool:
        return not any(self.c[1:])

    def __eq__(self, other) -> bool:
        if isinstance(other, FieldElement):
            return self.field == other.field and self.c == other.c
        if isinstance(other, (int, Fraction)):
            return self.is_rational() and self.c[0] == other
        return NotImplemented

    def __hash__(self) -> int:
        if self._hash is None:
            self._hash = hash(self.c)
        return self._hash

    def __float__(self) -> float:
        import math

        total = 0.0
        for m, a in enumerate(self.c):
            if a:
                total += float(a) * math.sqrt(self.field._prod[m])
        return total

    def __repr__(self) -> str:
        parts = []
        for m, a in enumerate(self.c):
            if a:
                parts.append(f"{a}" if m == 0 else f"{a}*{self.field._names[m]}")
        return " + ".join(parts) if parts else "0"


# The field that carries every construction in the literature we reproduce:
# sqrt(3) for the triangular (Eisenstein) lattice, sqrt(11) for the Moser
# spindle rotation arccos(5/6).  sqrt(33) comes along as their product.
QSQRT3_11 = Field((3, 11))


def embed(element: "FieldElement", target: Field) -> "FieldElement":
    """Re-express an element of one multiquadratic field inside a larger one.

    Needed because a spindle's rotation may require a square root the graph's
    own field does not have.  That is never a reason to skip the spindle: the
    field is a choice, and one more generator is free.  Whether a pair is
    forced monochromatic is a combinatorial fact about the graph that no field
    enters into -- the field only has to be big enough to write the rotated
    copy down afterwards.
    """
    src = element.field
    if src == target:
        return element
    coeffs = [0] * target.dim
    for m, c in enumerate(element.c):
        if not c:
            continue
        radicand = src._prod[m]
        placed = False
        for t in range(target.dim):
            if target._prod[t] == radicand:
                coeffs[t] = c
                placed = True
                break
        if not placed:
            raise ValueError(f"{target} does not contain sqrt({radicand})")
    return target.element(coeffs)
