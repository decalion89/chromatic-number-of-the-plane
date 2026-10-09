#!/usr/bin/env python3
"""Independent checks of the elementary number-theoretic claims introduced by the diff.
Standard library only.  Prints PASS/FAIL lines; exit status 1 on any failure."""
import sys
from fractions import Fraction as Fr

FAILS = []


def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""))
    if not cond:
        FAILS.append(name)


def squarefree(n):
    n = abs(n)
    d = 2
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 1
    return True


def is_square_Qp(u, p):
    """is the nonzero integer u a square in Q_p ? (p = 2 or odd p)"""
    v = 0
    while u % p == 0:
        u //= p
        v += 1
    if v % 2:
        return False
    if p == 2:
        return u % 8 == 1
    return pow(u % p, (p - 1) // 2, p) == 1


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def i_in_every_completion_above(m, p):
    """F = Q(sqrt m), m squarefree, not a square.  Every completion above p is Q_p (if m is a square in
    Q_p, split case) or Q_p(sqrt m).  i lies in it iff -1 or -m is a square in Q_p (resp. -1 in Q_p)."""
    if is_square_Qp(m, p):
        return is_square_Qp(-1, p)
    return is_square_Qp(-1, p) or is_square_Qp(-m, p)


def split_above_2_and_3(m):
    return i_in_every_completion_above(m, 2) and i_in_every_completion_above(m, 3)


# ---- Q(sqrt -73) and the real quadratic criterion
check("Q(sqrt -73): i lies in every completion above 2", i_in_every_completion_above(-73, 2))
check("Q(sqrt -73): i lies in every completion above 3", i_in_every_completion_above(-73, 3))
check("-73 is squarefree and -73 = 23 mod 24", squarefree(73) and (-73) % 24 == 23)
neg = [m for m in range(-1, -400, -1) if squarefree(m) and split_above_2_and_3(m)]
check("imaginary quadratic fields split above 2 and 3 (|m|<400): m = 23 mod 24; first after Q(i) is -73",
      all(m % 24 == 23 for m in neg) and neg[:2] == [-1, -73], str(neg))
pos_ok = all(split_above_2_and_3(d) == (d % 24 == 23) for d in range(2, 2000) if squarefree(d))
check("real quadratic Q(sqrt d), d squarefree < 2000: split above 2 and 3 iff d = 23 mod 24", pos_ok)

# ---- Q(sqrt2, sqrt7)
check("2 is a square mod 7 (7 splits in Q(sqrt2)); sqrt2 in Q_7", legendre(2, 7) == 1 and is_square_Qp(2, 7))
check("so Q(sqrt2,sqrt7) embeds in Q_7(sqrt7), residue field F_7 (place above 7 of residue degree 1)",
      is_square_Qp(2, 7) and not is_square_Qp(7, 7))
check("Prop. 2 hypotheses for Q(sqrt2,sqrt7): a=2 = 2 mod 3, b=7 = 7 mod 8", 2 % 3 == 2 and 7 % 8 == 7)
check("Prop. 3, p=7: 2a+3b on mu_8 = {+-1,+-i,+-2+-2i} takes values 2,5,3,4,3,5,2,4",
      [(2 * a + 3 * b) % 7 for (a, b) in [(1, 0), (-1, 0), (0, 1), (0, -1), (2, 2), (2, -2), (-2, 2), (-2, -2)]]
      == [2, 5, 3, 4, 3, 5, 2, 4]
      and all((a * a + b * b) % 7 == 1 for (a, b) in [(1, 0), (-1, 0), (0, 1), (0, -1), (2, 2), (2, -2), (-2, 2), (-2, -2)]))

# ---- finite witnesses for local fields
check("sqrt7 in Q_3 (7 = 1 mod 3)", is_square_Qp(7, 3))
check("sqrt11 in Q_7 (11 = 2^2 mod 7)", is_square_Qp(11, 7))
check("Q(sqrt7): (a) fails (7 = 3 mod 4), (b) holds (7 = 1 mod 3), so chi = 3", 7 % 4 == 3 and 7 % 3 != 2)
check("sqrt3 not in Q_3 (no unit triangle in Q_3^2)", not is_square_Qp(3, 3))

# ---- d = 47 and the list 47, 143, 167, 215
lst = [d for d in range(2, 216) if squarefree(d) and d % 24 == 23 and legendre(d, 7) == -1]
check("squarefree d<=215 with d=23 mod 24 and (d/7)=-1: 47,143,167,215", lst == [47, 143, 167, 215], str(lst))
check("examples: (23/7),(71/7) != -1; 59, 83 have (d/7)=-1, d=11 mod 24",
      legendre(23, 7) != -1 and legendre(71, 7) != -1 and all(legendre(d, 7) == -1 and d % 24 == 11 for d in (59, 83)))
check("47 = 3^2 mod 19 (sqrt47 in Q_19)", 47 % 19 == 9 and is_square_Qp(47, 19))
check("5 and 7 are squares mod 19 (Q(sqrt5,sqrt7) embeds in Q_19)", is_square_Qp(5, 19) and is_square_Qp(7, 19))
check("3 and -7 are not squares mod 19 (but -3 is)", legendre(3, 19) == -1 and legendre(-7, 19) == -1
      and legendre(-3, 19) == 1)
check("FC Cor. 2 open classes 47,143,167 mod 168 = {d=23 mod 24, (d/7)=-1}",
      sorted({d % 168 for d in range(1, 168 * 5) if d % 24 == 23 and legendre(d, 7) == -1}) == [47, 143, 167])


# ---- kappa_1 for f = 1 (p = 3 mod 4):  kappa_1 = max_{n in F_p^*} min_{s: n - s^2 in {0} u squares} ||2s/p||
def kappa1(p):
    sq = {(x * x) % p for x in range(p)}  # includes 0
    best = Fr(0)
    for n in range(1, p):
        m = None
        for u in range(0, p // 2 + 1):  # |2s| = u (centered), 2s = +-u
            for su in {u, (-u) % p}:
                s = (su * pow(2, p - 2, p)) % p
                if (n - s * s) % p in sq:
                    m = u
                    break
            if m is not None:
                break
        best = max(best, Fr(m, p))
    return best


def primes(lo, hi):
    out = []
    for n in range(lo, hi):
        if n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1)):
            out.append(n)
    return out


check("kappa_1(Q_3) = 1/3", kappa1(3) == Fr(1, 3))
check("kappa_1(Q_7) = 2/7", kappa1(7) == Fr(2, 7))
check("kappa_1(Q_19) = 4/19", kappa1(19) == Fr(4, 19))
vals = {p: kappa1(p) for p in primes(11, 3000) if p % 4 == 3}
mx = max(vals.values())
check("max kappa_1 over primes p = 3 mod 4, 11 <= p < 3000 is 4/19 (attained only at 19), all < 1/4",
      mx == Fr(4, 19) and [p for p in vals if vals[p] == mx] == [19] and all(v < Fr(1, 4) for v in vals.values()),
      "%d primes" % len(vals))
f = lambda p: Fr(p * p - 1, 2 * p)
import math
check("Prop. 6 bound at p=1001: (p^2-1)/(2p) - 2 sqrt p (ln p + 1) > 0",
      (1001 ** 2 - 1) / (2 * 1001) - 2 * math.sqrt(1001) * (math.log(1001) + 1) > 0,
      "%.4f" % ((1001 ** 2 - 1) / (2 * 1001) - 2 * math.sqrt(1001) * (math.log(1001) + 1)))


# ---- bonus: Lemma 19 for f = 3 directly: no F_49-linear map F_{7^6} -> F_49 maps mu_344 into A'_7
# A'_7 = {e in F_49 : e^8 = -1} (Lemma 13(i)), intrinsic, so no coordinates on F_49 are needed.
def lemma19_direct(f):
    p, deg = 7, 2 * f
    q2 = p ** deg
    order = q2 - 1

    def polymulmod(a, b, mod):
        res = [0] * (2 * deg - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    res[i + j] = (res[i + j] + x * y) % p
        for k in range(len(res) - 1, deg - 1, -1):
            c = res[k]
            if c:
                for j in range(deg + 1):
                    res[k - deg + j] = (res[k - deg + j] - c * mod[j]) % p
        return res[:deg]

    def pf(n):
        out, d = set(), 2
        while d * d <= n:
            while n % d == 0:
                out.add(d)
                n //= d
            d += 1
        if n > 1:
            out.add(n)
        return out

    def powx(e, mod):
        r = [1] + [0] * (deg - 1)
        b = [0, 1] + [0] * (deg - 2)
        while e:
            if e & 1:
                r = polymulmod(r, b, mod)
            b = polymulmod(b, b, mod)
            e >>= 1
        return r

    one = [1] + [0] * (deg - 1)
    # find a primitive polynomial (x has order q2-1)
    import itertools
    mod = None
    for coeffs in itertools.product(range(p), repeat=deg):
        cand = list(coeffs) + [1]  # monic, constant term first
        if cand[0] == 0:
            continue
        if powx(order, cand) != one:
            continue
        if all(powx(order // r, cand) != one for r in pf(order)):
            mod = cand
            break
    # exp / log tables
    enc = lambda v: sum(c * p ** i for i, c in enumerate(v))
    exp = [0] * order
    cur = one[:]
    xpoly = [0, 1] + [0] * (deg - 2)
    for k in range(order):
        exp[k] = enc(cur)
        cur = polymulmod(cur, xpoly, mod)
    log = {v: k for k, v in enumerate(exp)}
    assert len(log) == order

    def addenc(a, b):
        r, m = 0, 1
        while a or b:
            r += ((a % p + b % p) % p) * m
            a //= p
            b //= p
            m *= p
        return r

    minus_one_log = order // 2
    sub = p ** 2  # F_49 inside F_{7^deg}; Tr_{F_{q^2}/F_49}(w) = sum_{j<f} w^{49^j}
    flag = [False] * order
    for k in range(order):
        t = 0
        for j in range(f):
            t = addenc(t, exp[(k * sub ** j) % order])
        if t != 0 and (8 * log[t]) % order == minus_one_log:
            flag[k] = True
    qq = p ** f
    step = order // (qq + 1)  # mu_{q+1} = <g^step>
    good = [m for m in range(step) if all(flag[(m + step * l) % order] for l in range(qq + 1))]
    return good, step


good1, step1 = lemma19_direct(1)
check("sanity f=1: some F_49-linear map F_49 -> F_49 maps mu_8 into A'_7 (exactly one class of b mod mu_8)",
      len(good1) == 1, "classes of b: %d of %d" % (len(good1), step1))
good3, step3 = lemma19_direct(3)
check("Lemma 19, f=3, directly: no F_49-linear map F_{7^6} -> F_49 maps mu_344 into A'_7",
      len(good3) == 0, "%d classes of b tested (b up to mu_344), b=0 trivially fails" % step3)

print()
print("FAILED:", FAILS if FAILS else "none")
sys.exit(1 if FAILS else 0)
