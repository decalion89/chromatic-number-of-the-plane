"""b is determined by a: the chain is one squareness test, not a search.

Writing c = 1 + a, the closing condition |1 + a + b|^2 = 1/3 with |b| = 1 gives
c.bbar + cbar.b = -2/3 - |c|^2 =: T, and bbar = 1/b turns that into

    cbar . b^2  -  T . b  +  c  =  0,

a quadratic over K whose roots automatically have modulus one when its
discriminant is negative (their product is c/cbar).  So for EVERY unit step a
there are exactly two closing partners b in C, and the only question is
whether they lie in K:

    a extends to a chain  <=>  Delta = T^2 - 4|c|^2  is a square in K.

That replaces the pair sweep -- which found nothing, because independently
generated steps have wildly incompatible denominators -- with one test per a.

The test is done by reduction.  67, 199 and 331 are all 1 mod 33, so each
splits completely in Q(zeta_33) and reduction is just evaluation at a
primitive 33rd root of unity in F_l.  A square in K is a square in every
residue field, so a single non-residue rejects a outright.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.cyclotomic import CycloField

K = CycloField(33)
D = K.degree
ONE = K.rational(1)


def inv(a):
    rows = [list(K.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        s = Fraction(1) / M[c][c]
        M[c] = [v * s for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


Z11, Z3 = K.zeta(3), K.zeta(11)
Z6 = K.neg(K.mul(Z3, Z3))
g, p = K.zero(), ONE
for k in range(1, 11):
    p = K.mul(p, Z11) if k > 1 else Z11
    g = K.add(g, p if k in {1, 3, 4, 5, 9} else K.neg(p))
RHO = tuple(Fraction(1, 6) * x for x in K.add(K.rational(5), g))
assert K.mul(g, g) == K.rational(-11) and K.norm2(RHO) == ONE

t0 = time.time()
extra = set()
for coeffs in itertools.product(range(-1, 2), repeat=6):
    if not any(coeffs):
        continue
    a, q = K.zero(), ONE
    for c in coeffs:
        if c:
            a = K.add(a, tuple(Fraction(c) * x for x in q))
        q = K.mul(q, Z11)
    try:
        u = K.mul(a, inv(K.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if K.norm2(u) == ONE:
        extra.add(u)
zp, z = [], ONE
for _ in range(6):
    zp.append(z)
    z = K.mul(z, Z6)
steps = set()
for c in zp:
    for e in {ONE} | extra:
        base = K.mul(c, e)
        for m in (ONE, RHO, K.conj(RHO), K.mul(RHO, RHO),
                  K.mul(K.conj(RHO), K.conj(RHO))):
            steps.add(K.mul(base, m))
steps = sorted(steps)
print(f"{len(steps)} unit steps  [{time.time()-t0:.0f}s]", flush=True)

PRIMES = [67, 199, 331, 397, 463]
roots = {}
for L in PRIMES:
    for x in range(2, L):
        if pow(x, (L - 1) // 33, L) != 1:
            zt = pow(x, (L - 1) // 33, L)
            if pow(zt, 33, L) == 1 and all(pow(zt, 33 // q, L) != 1
                                           for q in (3, 11)):
                roots.setdefault(L, []).append(zt)
                if len(roots[L]) >= 3:
                    break
print("reduction primes: "
      + ", ".join(f"{L} ({len(roots.get(L, []))} roots)" for L in PRIMES),
      flush=True)


def is_square_mod(val):
    """False as soon as val is a non-residue in one residue field."""
    for L, zs in roots.items():
        for zt in zs:
            acc, ok = 0, True
            for i, x in enumerate(val):
                if x.denominator % L == 0:
                    ok = False
                    break
                acc = (acc + x.numerator * pow(x.denominator, -1, L)
                       * pow(zt, i, L)) % L
            if not ok:
                continue
            if acc == 0:
                continue
            if pow(acc, (L - 1) // 2, L) != 1:
                return False
    return True


survivors = []
TWO_THIRD = K.rational(Fraction(-2, 3))
for i, a in enumerate(steps):
    m = K.add(K.rational(2), K.add(a, K.conj(a)))          # |1 + a|^2
    T = K.sub(TWO_THIRD, m)
    delta = K.sub(K.mul(T, T), tuple(4 * x for x in m))
    if is_square_mod(delta):
        survivors.append((i, a, delta))
    if i % 1500 == 0:
        print(f"  ... {i}/{len(steps)}, {len(survivors)} survivors "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"{len(survivors)} of {len(steps)} steps survive the squareness filter "
      f"[{time.time()-t0:.0f}s]", flush=True)
for i, a, d in survivors[:6]:
    print(f"  survivor step {i}", flush=True)
