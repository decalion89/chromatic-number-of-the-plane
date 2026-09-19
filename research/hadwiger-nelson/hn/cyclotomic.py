"""Points in Z[zeta_n], where the rotations of odd order live.

`hn.transversal` ends on a fact about fields rather than about the plane: a
rotation of order n exists over a multiquadratic field exactly when n divides
24, the odd divisors of 24 are 1 and 3, and so the only odd cycle a
same-distance disjunction can ever be trapped in is the triangle -- which
traps two targets.  Every configuration in this package, de Grey's included,
lives over Q(sqrt3, sqrt5, sqrt7, sqrt11).  The ceiling came with the field.

So this module builds points somewhere else: the ring of integers Z[zeta_n] of
the cyclotomic field Q(zeta_n), whose Galois group is cyclic of order phi(n)
and which is multiquadratic only for n | 24.

Two facts make it a unit-distance world at all.  If z in Z[zeta_n] has |z| = 1
then z * conj(z) = 1 exactly, so every Galois conjugate has modulus 1 and
Kronecker's theorem makes z a root of unity.  The unit steps are therefore
precisely the 2n elements ±zeta_n^k -- a finite, exactly known step set, which
is all the search needs.  And Z[zeta_n] has rank phi(n) over Z, so for
phi(n) > 2 it is *dense* in the plane rather than discrete: the Eisenstein
lattice Z[omega] that every construction here has been built from is the
n = 3 case, rank 2, and the only one that is a lattice in the plane.

n = 15 is the first interesting choice.  Q(zeta_15) contains omega, so unit
triangles still exist -- the plane's clique number has to survive -- and it
contains zeta_5, so rotations of order 15 and 5 exist and the reachable odd
cycles become 3, 5 and 15, trapping up to 8 targets instead of 2.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import cos, gcd, pi, sin
from typing import List, Sequence, Tuple


def _poly_divmod(num: List[int], den: List[int]) -> List[int]:
    """Exact quotient of two integer polynomials, lowest degree first."""
    num = list(num)
    q = [0] * (len(num) - len(den) + 1)
    for i in range(len(q) - 1, -1, -1):
        c = num[i + len(den) - 1] // den[-1]
        q[i] = c
        if c:
            for j, d in enumerate(den):
                num[i + j] -= c * d
    return q


@lru_cache(maxsize=None)
def cyclotomic_poly(n: int) -> Tuple[int, ...]:
    """Phi_n, lowest degree first, from x^n - 1 divided by the Phi_d, d | n."""
    if n == 1:
        return (-1, 1)
    num = [0] * (n + 1)
    num[0], num[n] = -1, 1
    for d in range(1, n):
        if n % d == 0:
            num = _poly_divmod(num, list(cyclotomic_poly(d)))
    return tuple(num)


class CycloRing:
    """Z[zeta_n] in the power basis 1, zeta, ..., zeta^(phi(n)-1)."""

    def __init__(self, n: int):
        if n < 1:
            raise ValueError("n must be positive")
        self.n = n
        phi = cyclotomic_poly(n)
        self.degree = len(phi) - 1
        # x^degree == -(phi[0] + phi[1] x + ... ), phi monic
        self._fold = tuple(-c for c in phi[: self.degree])
        self._zeta_pow = tuple(self._reduce_monomial(k) for k in range(2 * n))

    def _reduce_monomial(self, k: int) -> Tuple[int, ...]:
        out = [0] * self.degree
        if k < self.degree:
            out[k] = 1
            return tuple(out)
        cur = [0] * self.degree
        cur[self.degree - 1] = 1
        for _ in range(k - self.degree + 1):
            top = cur[self.degree - 1]
            cur = [0] + cur[: self.degree - 1]
            if top:
                for j, f in enumerate(self._fold):
                    cur[j] += top * f
        return tuple(cur)

    # -- arithmetic ---------------------------------------------------------
    def zero(self) -> Tuple[int, ...]:
        return (0,) * self.degree

    def one(self) -> Tuple[int, ...]:
        return self._zeta_pow[0]

    def zeta(self, k: int = 1) -> Tuple[int, ...]:
        return self._zeta_pow[k % self.n]

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def mul(self, a, b):
        d = self.degree
        conv = [0] * (2 * d - 1)
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if y:
                    conv[i + j] += x * y
        out = list(conv[:d])
        for k in range(2 * d - 2, d - 1, -1):
            c = conv[k]
            if not c:
                continue
            for j, f in enumerate(self._zeta_pow[k]):
                out[j] += c * f
        return tuple(out)

    def conj(self, a):
        """Complex conjugation: zeta^j -> zeta^(-j)."""
        out = [0] * self.degree
        for j, x in enumerate(a):
            if not x:
                continue
            for i, f in enumerate(self._zeta_pow[(self.n - j) % self.n]):
                out[i] += x * f
        return tuple(out)

    # -- geometry -----------------------------------------------------------
    def norm2(self, a):
        """|a|^2 as a ring element.  Real, but kept in the power basis."""
        return self.mul(a, self.conj(a))

    def is_unit_apart(self, a, b) -> bool:
        return self.norm2(self.sub(a, b)) == self.one()

    def unit_steps(self) -> List[Tuple[int, ...]]:
        """The 2n elements of modulus 1: every unit step there is.

        |z| = 1 forces z * conj(z) = 1 as an algebraic identity, so all Galois
        conjugates have modulus 1 and Kronecker makes z a root of unity.  For
        odd n the 2n-th roots are already ±zeta_n^k; for even n, -1 is itself
        a power and the set closes at n.
        """
        seen = {}
        for k in range(self.n):
            z = self.zeta(k)
            seen[z] = None
            seen[self.neg(z)] = None
        return [z for z in seen if self.norm2(z) == self.one()]

    def to_complex(self, a) -> complex:
        t = 2 * pi / self.n
        return sum(x * complex(cos(t * j), sin(t * j)) for j, x in enumerate(a))


def rotation_orders(n: int) -> List[int]:
    """Rotation orders available over Q(zeta_n): the divisors of n, plus 2n
    when n is odd, since -zeta_n^k is then a primitive 2n-th root."""
    m = n if n % 2 == 0 else 2 * n
    return [d for d in range(1, m + 1) if m % d == 0]


def odd_cycle_reach(n: int) -> int:
    """Largest same-distance target set trappable over Q(zeta_n).

    Reads `hn.transversal` with this field's rotation orders instead of the
    multiquadratic ones.
    """
    from hn.transversal import same_distance_ceiling

    return same_distance_ceiling(rotation_orders(n))


class CycloField(CycloRing):
    """Q(zeta_n): the ring with rational coefficients.

    Points of a construction are not algebraic integers -- the Moser spindle's
    rotation is (5 + sqrt(-11))/6 -- so the coefficients have to be fractions.
    The arithmetic is otherwise identical, inherited unchanged.
    """

    def zero(self):
        return (Fraction(0),) * self.degree

    def one(self):
        return tuple(Fraction(c) for c in self._zeta_pow[0])

    def zeta(self, k: int = 1):
        return tuple(Fraction(c) for c in self._zeta_pow[k % self.n])

    def rational(self, q) -> Tuple:
        out = [Fraction(0)] * self.degree
        out[0] = Fraction(q)
        return tuple(out)

    def scale(self, a, q):
        q = Fraction(q)
        return tuple(x * q for x in a)

    def sqrt_disc(self, p: int):
        """The quadratic Gauss sum: sqrt(p*) for p* = p if p%4==1 else -p.

        g = sum_k (k|p) zeta_p^k equals sqrt(p) when p = 1 mod 4 and
        sqrt(-p) otherwise, which is exactly the radical a spindle rotation
        needs.  Requires p | n so that zeta_p is available here.
        """
        if self.n % p:
            raise ValueError(f"zeta_{p} is not in Q(zeta_{self.n})")
        step = self.n // p
        residues = {(k * k) % p for k in range(1, p)}
        out = [Fraction(0)] * self.degree
        for k in range(1, p):
            sign = 1 if k % p in residues else -1
            for j, c in enumerate(self._zeta_pow[(k * step) % self.n]):
                out[j] += sign * c
        return tuple(out)


def moser_rotation(field: "CycloField"):
    """Multiplication by (5 + sqrt(-11))/6: the Moser spindle's rotation.

    The spindle joins two points at distance sqrt(3) from a pivot, so its angle
    has cos = 1 - 1/6 = 5/6 and sin = sqrt(11)/6.  As a complex number that is
    (5 + i sqrt 11)/6 = (5 + sqrt(-11))/6, of modulus 1, and sqrt(-11) is the
    Gauss sum over zeta_11.  So the rotation is available in Q(zeta_n) for any
    n divisible by 11 -- no square root outside the field, and no separate
    real field for the coordinates.
    """
    r = field.sqrt_disc(11)
    five = field.rational(5)
    return field.scale(field.add(five, r), Fraction(1, 6))
