"""Choose the field instead of searching it: m of degree three.

The conic y1^2 = 3(4 - r^2) has the rational point (r, y1) = (1, 3), so it is
rationally parametrised by r = (m^2 - 6m - 3)/(m^2 + 3), and the second
condition becomes a single quartic,

    w^2 = -39 m^4 + 540 m^3 - 666 m^2 - 324 m + 297  =  3 N(m).

Two conditions come free.  r + 2 = 3(m-1)^2/(m^2+3) and 2 - r = (m+3)^2/(m^2+3)
are both non-negative, so |r| <= 2 at EVERY real embedding and the unit step a
always exists.  And a closing chain forces some prime above 3 to be moved by
conjugation, proved here, so constructing the chain settles the condition at 3
rather than needing it imposed.

What is left is the condition that decides blocking: every prime of F above 5
must have residue degree at least 3.  r inherits m's degree, so m is taken
cubic and F = Q(m, w) has degree 6 with defining polynomial
g(X) = prod_i (X^2 - 3N(m_i)).  All residue degrees are at least 3 exactly
when g mod 5 has no factor of degree 1 or 2, i.e. when gcd(g, X^25 - X) is
constant -- one polynomial gcd per candidate.
"""
import time


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def pmod(a, b, p):
    a = [x % p for x in a]
    b = trim([x % p for x in b])
    inv = pow(b[-1], -1, p)
    trim(a)
    while len(a) >= len(b):
        f = a[-1] * inv % p
        sh = len(a) - len(b)
        for i, y in enumerate(b):
            a[i + sh] = (a[i + sh] - f * y) % p
        trim(a)
    return a


def pmul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return trim(r)


def pgcd(a, b, p):
    a, b = trim([x % p for x in a]), trim([x % p for x in b])
    while b:
        a, b = b, pmod(a[:], b, p)
    return a


def xpow(e, g, p):
    """X^e mod g."""
    res, base = [1], pmod([0, 1], g, p)
    while e:
        if e & 1:
            res = pmod(pmul(res, base, p), g, p)
        base = pmod(pmul(base, base, p), g, p)
        e >>= 1
    return res


def min_factor_degree_at_least_3(g, p=5):
    """True when g mod p is squarefree with every factor of degree >= 3."""
    gm = trim([c % p for c in g])
    if len(gm) != len(g):
        return False
    d = [(i * c) % p for i, c in enumerate(gm)][1:]
    if len(pgcd(gm, d, p)) != 1:
        return False                      # not squarefree: 5 ramifies
    for k in (1, 2):
        xp = xpow(p ** k, gm, p)
        while len(xp) < 2:
            xp.append(0)
        diff = xp[:]
        diff[1] = (diff[1] - 1) % p
        if len(pgcd(gm, diff, p)) > 1:
            return False
    return True


def real_cubic_roots(c2, c1, c0):
    import math
    p = c1 - c2 * c2 / 3.0
    q = 2 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
    if p >= 0:
        return None
    disc = -(4 * p ** 3 + 27 * q * q)
    if disc <= 1e-9:
        return None
    arg = 3 * q / (2 * p) * math.sqrt(-3.0 / p)
    if abs(arg) > 1:
        return None
    return sorted(2 * math.sqrt(-p / 3.0)
                  * math.cos(math.acos(arg) / 3.0 - 2 * math.pi * k / 3.0)
                  - c2 / 3.0 for k in range(3))


def has_rational_root(c2, c1, c0):
    for d in range(1, abs(c0) + 2):
        if c0 % d == 0:
            for s in (d, -d):
                if s ** 3 + c2 * s * s + c1 * s + c0 == 0:
                    return True
    return c0 == 0


def N3(m):
    """V(m): the quartic itself, which already carries the factor 3.

    An earlier pass multiplied it by 3 again.  Q(sqrt(3V)) is not Q(sqrt V),
    so that searched the wrong fields entirely -- caught by writing out the
    derivation numerically, where |a| = |b| = 1 and |1 + a + b|^2 = 1/3 hold
    to machine precision for every real m with V(m) > 0.
    """
    return -39 * m ** 4 + 540 * m ** 3 - 666 * m ** 2 - 324 * m + 297


t0 = time.time()
good, tried = [], 0
for c2 in range(-9, 10):
    for c1 in range(-24, 25):
        for c0 in range(-24, 25):
            rs = real_cubic_roots(c2, c1, c0)
            if rs is None or has_rational_root(c2, c1, c0):
                continue
            vals = [N3(m) for m in rs]
            if min(vals) <= 1e-9:
                continue                   # F must be totally real
            tried += 1
            g = [1.0]
            for v in vals:
                g = [0.0, 0.0] + g
                for i in range(len(g) - 2):
                    g[i] -= v * g[i + 2]
            gi = []
            ok = True
            for x in g:
                n = round(x)
                if abs(x - n) > 1e-3 * max(1.0, abs(x)):
                    ok = False
                    break
                gi.append(n)
            if not ok or len(gi) != 7:
                continue
            if min_factor_degree_at_least_3(gi):
                good.append((c2, c1, c0, gi))
                print(f"  *** m root of x^3 + {c2}x^2 + {c1}x + {c0}: every "
                      f"prime above 5 has residue degree >= 3", flush=True)
                if len(good) >= 40:
                    break
        if len(good) >= 40:
            break
    if len(good) >= 40:
        break
print(f"{tried} totally real cubics with 3N(m) totally positive; "
      f"{len(good)} give a blocking-capable field  [{time.time()-t0:.0f}s]",
      flush=True)
