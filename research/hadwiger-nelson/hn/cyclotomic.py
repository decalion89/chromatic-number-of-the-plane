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

    def _modulus(self):
        """Phi_n as rational coefficients, for polynomial inversion."""
        return [Fraction(c) for c in cyclotomic_poly(self.n)]

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


def _poly_inv_mod(a: List[Fraction], mod: List[Fraction]) -> List[Fraction]:
    """Inverse of a modulo mod, over Q, by the extended Euclidean algorithm."""
    def deg(p):
        d = len(p) - 1
        while d >= 0 and p[d] == 0:
            d -= 1
        return d

    def sub_scaled(p, q, c, shift):
        out = list(p)
        for i, x in enumerate(q):
            while len(out) <= i + shift:
                out.append(Fraction(0))
            out[i + shift] -= c * x
        return out

    r0, r1 = list(mod), list(a)
    s0, s1 = [Fraction(0)], [Fraction(1)]
    while deg(r1) >= 0:
        d0, d1 = deg(r0), deg(r1)
        if d0 < d1:
            r0, r1, s0, s1 = r1, r0, s1, s0
            continue
        c = r0[d0] / r1[d1]
        r0 = sub_scaled(r0, r1, c, d0 - d1)
        s0 = sub_scaled(s0, s1, c, d0 - d1)
        if deg(r0) < deg(r1):
            r0, r1, s0, s1 = r1, r0, s1, s0
    if deg(r0) != 0:
        raise ZeroDivisionError("element is not invertible")
    scale = r0[0]
    return [x / scale for x in s0]


def unit_polygon(field: "CycloField", order: int) -> List[Tuple]:
    """The regular `order`-gon of side 1, centred at the origin.

    Its circumradius is the magic radius 1/(2 sin(pi/order)), because a chord
    of a circle subtending one full rotation step has length 1 exactly there.
    And it has a closed form needing no radicals at all: with z the primitive
    root of that order,

        v_k = z^k / (z - 1),

    since |v_1 - v_0| = |z - 1| / |z - 1| = 1.  Adjacent vertices are one
    apart, so the polygon *is* a unit-distance cycle of length `order` -- odd
    when the order is odd, which is the whole point.
    """
    if field.n % order:
        raise ValueError(f"zeta_{order} is not in Q(zeta_{field.n})")
    z = field.zeta(field.n // order)
    denom = field.sub(z, field.one())
    inv = tuple(_poly_inv_mod(list(denom), list(field._modulus()))[: field.degree]
                + [Fraction(0)] * field.degree)[: field.degree]
    return [field.mul(field.zeta((field.n // order) * k), inv) for k in range(order)]


def unit_steps_extended(field: "CycloField", rot, m_max: int = 2) -> List[Tuple]:
    """Unit steps beyond the roots of unity: ±zeta^k * rot^m.

    `unit_steps` returns the algebraic *integers* of modulus 1, which Kronecker
    pins to the 2n roots of unity.  But a unit step only has to have modulus 1,
    and Q(zeta_n) has far more such elements than Z[zeta_n] does -- the Moser
    rotation (5 + sqrt(-11))/6 is one, with a denominator of 6.

    Multiplying the roots of unity by powers of such a rotation gives a much
    larger step set, and a richer one: a step set closed under multiplication
    by zeta_n forces every ball grown from it to be rotation-invariant, and
    invariance is what stops a core shrinking below a full orbit.  Powers of
    rot break that closure while keeping every step exactly unit length.

    It also puts spindles everywhere rather than only about the origin, which
    is what a rotation applied as a generator, rather than as a step, can never
    do.
    """
    if field.norm2(rot) != field.one():
        raise ValueError("rot must have modulus 1")
    powers = [field.one()]
    for _ in range(m_max):
        powers.append(field.mul(powers[-1], rot))
    inv = field.conj(rot)                       # modulus 1, so the inverse
    back = [field.one()]
    for _ in range(m_max):
        back.append(field.mul(back[-1], inv))
    seen = {}
    for m in powers + back[1:]:
        for k in range(field.n):
            z = field.mul(field.zeta(k), m)
            seen[z] = None
            seen[field.neg(z)] = None
    return [z for z in seen if field.norm2(z) == field.one()]


def unit_steps_from_quotients(field: "CycloField", spread: int = 1,
                              max_den: int = 200, limit: int = 400) -> List[Tuple]:
    """Modulus-one elements built as alpha / conj(alpha).

    `unit_steps` gives the algebraic integers of modulus one, which Kronecker
    pins to the 2n roots of unity, and those are closed under multiplying by
    zeta_n -- so a ball grown from them is rotation-invariant, and in an
    invariant ball no proper subset of a target orbit is ever forced.

    Hilbert 90 says every modulus-one element of Q(zeta_n) is alpha/conj(alpha)
    for some alpha, and those are *not* closed under zeta_n.  They are also
    where the spindles come from: the Moser rotation (5 + sqrt(-11))/6 is one
    of them, and a field without it needs its own.  Q(zeta_15) has no sqrt(-11)
    and its integer walk is 3-chromatic, which is the whole problem.

    `spread` bounds the coefficients of alpha, `max_den` the denominator of the
    quotient, since an unbounded one makes the integer walk overflow rather
    than help.
    """
    from itertools import product

    out = {}
    rng = range(-spread, spread + 1)
    for coeffs in product(rng, repeat=min(4, field.degree)):
        if not any(coeffs):
            continue
        a = list(field.zero())
        for j, c in enumerate(coeffs):
            a[j] = Fraction(c)
        a = tuple(a)
        try:
            inv = _poly_inv_mod(list(field.conj(a)), field._modulus())
        except ZeroDivisionError:
            continue
        inv = tuple((inv + [Fraction(0)] * field.degree)[: field.degree])
        q = field.mul(a, inv)
        if field.norm2(q) != field.one():
            continue
        den = 1
        for x in q:
            den = den * x.denominator // __import__("math").gcd(den, x.denominator)
        if den > max_den:
            continue
        out[q] = None
        out[field.neg(q)] = None
        if len(out) >= limit:
            break
    return list(out)
