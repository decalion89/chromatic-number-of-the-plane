"""Q(zeta_21, sqrt-11): the smallest field with triangles, blocking and a spindle.

The theorem said it. Triangles need zeta_6, so 3 | n. Blocking needs a prime
splitting completely in Q(zeta_7), so 7 | n. And the Moser spindle needs
sqrt(-11), which no Q(zeta_21) rotation supplies, because the Eisenstein
spindle condition 12m - 1 = 3t^2 or 7t^2 is impossible modulo 12.

Q(zeta_21) has degree 12 and adjoining sqrt(-11) doubles it. Elements are
a + b s with a, b in Q(zeta_21) and s^2 = -11, so twenty-four rational
coordinates. Complex conjugation fixes the rationals, sends zeta to its
inverse, and negates s since s is purely imaginary:

    conj(a + b s) = conj(a) - conj(b) s.

That is all the arithmetic a unit-distance graph needs: |z|^2 = z conj(z),
and "one apart" is |z - w|^2 == 1 exactly.
"""
import sys
from fractions import Fraction
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField

K = CycloField(21)
D = K.degree


class Ext:
    """a + b*sqrt(-11) over Q(zeta_21)."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=None):
        self.a = a
        self.b = b if b is not None else K.zero()

    def __eq__(self, o):
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return f"Ext({self.a!r}, {self.b!r})"


def e_zero():
    return Ext(K.zero(), K.zero())


def e_of(a):
    return Ext(a, K.zero())


def e_add(x, y):
    return Ext(K.add(x.a, y.a), K.add(x.b, y.b))


def e_sub(x, y):
    return Ext(K.sub(x.a, y.a), K.sub(x.b, y.b))


def e_neg(x):
    return Ext(K.neg(x.a), K.neg(x.b))


def e_mul(x, y):
    # (a + bs)(c + ds) = (ac - 11bd) + (ad + bc)s
    ac = K.mul(x.a, y.a)
    bd = K.mul(x.b, y.b)
    return Ext(K.sub(ac, tuple(Fraction(11) * t for t in bd)),
               K.add(K.mul(x.a, y.b), K.mul(x.b, y.a)))


def e_conj(x):
    return Ext(K.conj(x.a), K.neg(K.conj(x.b)))


def e_norm2(x):
    """|z|^2 as a rational, or None if it is not rational."""
    z = e_mul(x, e_conj(x))
    if any(t for t in z.b):
        return None
    c = z.a
    if any(t for t in c[1:]):
        return None
    return c[0]


ONE = e_of(K.rational(1))
S11 = Ext(K.zero(), K.rational(1))
assert e_norm2(S11) == 11, "sqrt(-11) has modulus squared 11"
Z3 = e_of(K.zeta(7))
Z6 = e_neg(e_mul(Z3, Z3))
assert e_norm2(Z6) == 1 and e_norm2(e_sub(Z6, ONE)) == 1, "zeta_6 is a sixth root"
Z7 = e_of(K.zeta(3))
assert e_norm2(Z7) == 1

# the Moser rotation (5 + sqrt(-11)) / 6
SIX = e_of(K.rational(Fraction(1, 6)))
RHO = e_mul(SIX, e_add(e_of(K.rational(5)), S11))
print(f"degree 24 field built; |rho|^2 = {e_norm2(RHO)}, "
      f"|1 - rho|^2 = {e_norm2(e_sub(ONE, RHO))}", flush=True)
assert e_norm2(RHO) == 1
assert e_norm2(e_sub(ONE, RHO)) == Fraction(1, 3), "the Moser chord"

# the spindle: two unit rhombi sharing the origin
rh = [e_zero(), ONE, Z6, e_add(ONE, Z6)]
pts = list(rh) + [e_mul(RHO, p) for p in rh[1:]]
seen, P = set(), []
for p in pts:
    if p not in seen:
        seen.add(p)
        P.append(p)
E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
     if e_norm2(e_sub(P[j], P[i])) == 1]
print(f"spindle over Q(zeta_21, sqrt-11): {len(P)} points, {len(E)} edges",
      flush=True)

from pysat.solvers import Solver
for k in (3, 4):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"  {k}-colourable: {r}"
          + ("" if r else "   *** chi >= 4 IN A BLOCKING FIELD ***"),
          flush=True)
    if r:
        break
