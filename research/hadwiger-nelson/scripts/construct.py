"""Construct the resonance instead of searching for it.

A union of a lattice with its rotated copy gains an edge exactly when
|rho p - q| = 1 for lattice points p, q.  Expanding, with z = p conj(q) and
N for the norm,

    Re(rho z) = (N(p) + N(q) - 1) / 2 =: R,

and since |rho| = 1 this is a quadratic: rho^2 z - 2R rho + conj(z) = 0, so

    rho = (R +- sqrt(R^2 - N(z))) / z.

R^2 - N(z) is usually negative, so rho lies in the field exactly when
sqrt(-(N(z) - R^2)) does.  Q(zeta_21) contains sqrt(-3) and sqrt(-7), so the
condition is that N(z) - R^2 be three or seven times a rational square.

That is a finite check over pairs of lattice points, and it either produces a
rotation that is not a root of unity -- a Moser analogue in a field that
blocks -- or shows that every solution is trivial.

Worth noting what the Moser rotation itself is in this language: p at sqrt 3,
q at 1, R = 3/2, N(z) = 3, and 3 - 9/4 = 3/4 is three times a square. The two
roots come out 1 and zeta_6^-1 -- both trivial. So the classical case needs
larger p and q, and so does anything new.
"""
import sys, time, itertools
from fractions import Fraction
from math import isqrt
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField
from pysat.solvers import Solver

F = CycloField(21)
D = F.degree
one = F.rational(1)
z3 = F.zeta(7)
z6 = F.neg(F.mul(z3, z3))
S3 = F.sub(F.add(z3, z3), F.neg(one))      # 2 zeta_3 + 1 = sqrt(-3)
assert F.mul(S3, S3) == F.rational(-3), "sqrt(-3) in the field"
z7 = F.zeta(3)
S7 = None
for e in range(1, 21):                      # sqrt(-7) from the zeta_7 part
    g = F.zero()
    for j in range(1, 7):
        t = F.zeta(3 * (j * j % 7))
        g = F.add(g, t)
    break
# Gauss sum: sum_{j} legendre(j,7) zeta_7^j = sqrt(-7)
S7 = F.zero()
for j in range(1, 7):
    leg = 1 if pow(j, 3, 7) == 1 else -1
    t = F.zeta(3 * j)
    S7 = F.add(S7, t if leg > 0 else F.neg(t))
assert F.mul(S7, S7) == F.rational(-7), "sqrt(-7) in the field"
print("sqrt(-3) and sqrt(-7) both confirmed in Q(zeta_21)", flush=True)


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        pr = next(r for r in range(c, D) if M[r][c])
        M[c], M[pr] = M[pr], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def eis(a, b):
    v = F.zero()
    for _ in range(abs(a)):
        v = F.add(v, one if a > 0 else F.neg(one))
    for _ in range(abs(b)):
        v = F.add(v, z6 if b > 0 else F.neg(z6))
    return v


R_ = 3
patch = sorted({eis(a, b) for a in range(-R_, R_ + 1)
                for b in range(-R_, R_ + 1)})
print(f"Eisenstein patch: {len(patch)} points", flush=True)

roots = set()
p_ = one
for _ in range(21):
    p_ = F.mul(p_, F.zeta(1))
    roots.add(p_)
    roots.add(F.neg(p_))


def sq_times(m, d):
    """Is m = d * t^2 for rational t?"""
    num, den = m.numerator * d, m.denominator * d
    q = Fraction(m, d)
    n, dd = q.numerator, q.denominator
    return isqrt(n) ** 2 == n and isqrt(dd) ** 2 == dd and n >= 0


found, tested, t0 = [], 0, time.time()
for p in patch:
    Np = F.norm2(p)
    if Np == F.zero():
        continue
    for q in patch:
        Nq = F.norm2(q)
        if Nq == F.zero():
            continue
        try:
            npv, nqv = Fraction(str(Np.c[0] if hasattr(Np, 'c') else Np)), None
        except Exception:
            pass
        npf = Fraction(Np[0]) if isinstance(Np, tuple) else None
        nqf = Fraction(Nq[0]) if isinstance(Nq, tuple) else None
        if npf is None or nqf is None:
            continue
        Rv = (npf + nqf - 1) / 2
        z = F.mul(p, F.conj(q))
        Nz = Fraction(F.norm2(z)[0])
        m = Nz - Rv * Rv
        if m <= 0:
            continue
        tested += 1
        for d, rad in ((3, S3), (7, S7)):
            if not sq_times(m, d):
                continue
            t = Fraction(m, d)
            tt = Fraction(isqrt(t.numerator), isqrt(t.denominator))
            for sign in (1, -1):
                num = F.add(F.rational(Rv),
                            tuple(Fraction(sign) * tt * x for x in rad))
                rho = F.mul(num, inverse(z))
                if F.norm2(rho) != one or rho in roots:
                    continue
                found.append(rho)
print(f"{tested} pairs with a positive discriminant; "
      f"{len(found)} non-trivial rotations  [{time.time()-t0:.0f}s]",
      flush=True)

# Rank by how many cross edges each rotation actually creates: two or four is
# not enough to force a fourth colour, and the Moser configuration needs a
# particular one.  Then try unions of the patch with SEVERAL rotated copies,
# since each adds its own cross edges to the same lattice.
ranked, seen = [], set()
for rho in found:
    if rho in seen:
        continue
    seen.add(rho)
    img = {F.mul(rho, x) for x in patch}
    bs = set(patch)
    cross = sum(1 for a in bs - img for b in img - bs
                if F.norm2(F.sub(a, b)) == one)
    ranked.append((cross, rho))
ranked.sort(key=lambda t: -t[0])
print(f"{len(ranked)} distinct rotations; cross edges "
      f"{[c for c, _ in ranked[:8]]}", flush=True)


def chi_of(P, hi=5):
    E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
         if F.norm2(F.sub(P[j], P[i])) == one]
    for kk in range(3, hi + 1):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(P))]
        for a, b in E:
            for c in range(kk):
                cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return kk, len(E)
    return hi + 1, len(E)


print("\nsingle copies, richest first:", flush=True)
for cross, rho in ranked[:6]:
    P = sorted(set(patch) | {F.mul(rho, x) for x in patch})
    k, e = chi_of(P)
    print(f"  {cross} cross edges: {len(P)} points, {e} edges, chi = {k}"
          + ("   *** 4-CHROMATIC ***" if k >= 4 else ""), flush=True)
    if k >= 4:
        raise SystemExit

print("\nseveral copies at once:", flush=True)
cur = set(patch)
for idx, (cross, rho) in enumerate(ranked[:14], 1):
    cur |= {F.mul(rho, x) for x in patch}
    if len(cur) > 700:
        break
    k, e = chi_of(sorted(cur))
    print(f"  {idx} rotations: {len(cur)} points, {e} edges, chi = {k}"
          + ("   *** 4-CHROMATIC OVER A BLOCKING FIELD ***" if k >= 4 else ""),
          flush=True)
    if k >= 4:
        break
