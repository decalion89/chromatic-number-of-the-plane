"""The chain of three is a genus-one curve, and it depends only on Re(a).

With c = 1 + a and m = |c|^2 = 2 + r, r = a + abar, the closing quadratic
cbar b^2 - T b + c = 0 has T = -2/3 - m and discriminant

    Delta = T^2 - 4m = r^2 + (4/3) r - 8/9,

which involves ONLY the real part of a.  K = F(sqrt-3), so Delta -- negative
at every real embedding -- is a square in K exactly when -Delta/3 is a square
in F, i.e. when 3(8 - 12r - 9r^2) is.  And a itself exists as a unit step with
a + abar = r exactly when 3(4 - r^2) is a square in F.  So the whole question
is

    y1^2 = 3(4 - r^2),     y2^2 = 3(8 - 12 r - 9 r^2),      r, y1, y2 in F:

two quadrics in P^3, a curve of genus one, whose F-points are the chains.

Rational r is no use for the wider programme -- it puts a and b in quadratic
fields, so the chain field is multiquadratic and can never block -- but it
tests the derivation, and the survivors of the reduction filter test it too.
"""
import sys, itertools, time, cmath
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField
from pysat.solvers import Solver

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
ONE_PLUS = K.add(ONE, Z6)


def val(t):
    return sum(complex(float(x)) * cmath.exp(2j * cmath.pi * i / 33)
               for i, x in enumerate(t))


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
nums = [val(s) for s in steps]
print(f"{len(steps)} unit steps  [{time.time()-t0:.0f}s]", flush=True)

PRIMES = [67, 199, 331, 397, 463]
roots = {}
for L in PRIMES:
    for x in range(2, L):
        zt = pow(x, (L - 1) // 33, L)
        if pow(zt, 33, L) == 1 and all(pow(zt, 33 // q, L) != 1 for q in (3, 11)):
            roots.setdefault(L, []).append(zt)
            if len(roots[L]) >= 3:
                break


def square_mod(v):
    for L, zs in roots.items():
        for zt in zs:
            acc, ok = 0, True
            for i, x in enumerate(v):
                if x.denominator % L == 0:
                    ok = False
                    break
                acc = (acc + x.numerator * pow(x.denominator, -1, L)
                       * pow(zt, i, L)) % L
            if ok and acc and pow(acc, (L - 1) // 2, L) != 1:
                return False
    return True


THREE = K.rational(3)
FOUR3 = K.rational(Fraction(4, 3))
EIGHT9 = K.rational(Fraction(8, 9))
PAIRS = [(46, 3929), (47, 6036)]


def rank_of(vs):
    rows = [[Fraction(c) for c in v] for v in vs]
    r = 0
    for c in range(D):
        pr = next((t for t in range(r, len(rows)) if rows[t][c]), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        f = rows[r][c]
        rows[r] = [x / f for x in rows[r]]
        for t in range(len(rows)):
            if t != r and rows[t][c]:
                k2 = rows[t][c]
                rows[t] = [x - k2 * y for x, y in zip(rows[t], rows[r])]
        r += 1
    return r


BASIS = [ONE, Z3, g, K.mul(Z3, g)]
BR = rank_of(BASIS)
from hn.homcol import has_homomorphism

for ia, ib in PAIRS:
    a, b = steps[ia], steps[ib]
    s3 = K.add(K.add(ONE, a), b)
    assert K.norm2(s3) == K.rational(Fraction(1, 3)), (ia, ib)
    out_a = rank_of(BASIS + [a]) > BR
    out_b = rank_of(BASIS + [b]) > BR
    ws = [ONE, a, b]
    pts, B = [K.zero()], K.zero()
    for w in ws:
        for off in (w, K.mul(Z6, w), K.mul(ONE_PLUS, w)):
            pts.append(K.add(B, off))
        B = K.add(B, K.mul(ONE_PLUS, w))
    uniq, seen = [], set()
    for q in pts:
        if q not in seen:
            seen.add(q); uniq.append(q)
    zs = [val(q) for q in uniq]
    E = [(i, j) for i in range(len(uniq)) for j in range(i + 1, len(uniq))
         if abs(abs(zs[i] - zs[j]) - 1) < 1e-7
         and K.norm2(K.sub(uniq[j], uniq[i])) == ONE]
    chi = None
    for kk in (3, 4, 5):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(uniq))]
        for x, y in E:
            for c in range(kk):
                cls.append([-(1 + x * kk + c), -(1 + y * kk + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            if sv.solve():
                chi = kk
                break
    vecs = []
    for x, y in E:
        vecs.append(K.sub(uniq[y], uniq[x]))
        vecs.append(K.sub(uniq[x], uniq[y]))
    dn = 1
    for v in vecs:
        for x in v:
            dn = dn * x.denominator // gcd(dn, x.denominator)
    iv = set()
    for v in vecs:
        w2 = tuple(int(x * dn) for x in v)
        gg = 0
        for t in w2:
            gg = gcd(gg, abs(t))
        iv.add(tuple(t // gg for t in w2) if gg > 1 else w2)
    rk = rank_of([tuple(Fraction(x) for x in v) for v in iv])
    phi, _ = has_homomorphism(sorted(iv), 5)
    print(f"chain ({ia},{ib}): {len(uniq)} points, {len(E)} edges, chi = {chi}, "
          f"{len(iv)} directions, module rank {rk}; a outside subfield "
          f"{out_a}, b outside {out_b}; "
          + ("BLOCKS" if phi is None else "coset colouring exists"), flush=True)
