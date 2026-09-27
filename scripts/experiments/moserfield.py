"""A field with de Grey's spindle AND the residue degree that blocks.

A distance d spindles iff K holds a rotation with |1 - rho|^2 = 1/d^2.  With
t = rho + rhobar = 2 - 1/d^2 one gets 4 - t^2 = (4d^2 - 1)/d^4, so since
K = F(sqrt-3) the condition is

    3(4d^2 - 1)  is a square in F.

At d^2 = 3 -- de Grey's distance, the rhombus tip -- that is 33.  And
F = Q(m)(sqrt V) contains Q(sqrt33) exactly when V = 33 . (a square in Q(m)),
by the usual argument: if sqrt V lies in Q(m, sqrt D) then V = a^2 + b^2 D +
2ab sqrt D forces ab = 0, and b = 0 would make V a square outright.

Crucially 33 is POSITIVE, so V = 33 s^2 is compatible with V totally positive,
which the chain needs.  Such a cubic would give

    F = Q(m, sqrt33),   K = Q(m, sqrt33, sqrt-3) = Q(m, sqrt-3, sqrt-11),

holding the Moser rotation (5 + sqrt-11)/6 -- so de Grey's route to five
colours is available -- while still being a degree-12 field whose residue
degrees at 5 can be 3 or more, which a multiquadratic field never is.

Searched here over totally real cubics, filtered first by the cheap necessary
condition that N(V)/33^3 be a rational square.
"""
import math, time
from fractions import Fraction as Fr

Q = (0, 0, 0)


def real_cubic_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
    if p >= 0:
        return None
    if -(4 * p ** 3 + 27 * q * q) <= 1e-9:
        return None
    arg = 3 * q / (2 * p) * math.sqrt(-3.0 / p)
    if abs(arg) > 1:
        return None
    return [2 * math.sqrt(-p / 3.0)
            * math.cos(math.acos(arg) / 3.0 - 2 * math.pi * k / 3.0) - c2 / 3.0
            for k in range(3)]


def has_rational_root(c2, c1, c0):
    if c0 == 0:
        return True
    for d in range(1, abs(c0) + 1):
        if c0 % d == 0:
            for s in (d, -d):
                if s ** 3 + c2 * s * s + c1 * s + c0 == 0:
                    return True
    return False


def V(m):
    return -39 * m ** 4 + 540 * m ** 3 - 666 * m ** 2 - 324 * m + 297


def issq(n):
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


t0 = time.time()
hits, tried = [], 0
for c2 in range(-14, 15):
    for c1 in range(-40, 41):
        for c0 in range(-40, 41):
            rs = real_cubic_roots(c2, c1, c0)
            if rs is None or has_rational_root(c2, c1, c0):
                continue
            vs = [V(m) for m in rs]
            if min(vs) <= 1e-9:
                continue
            tried += 1
            nv = vs[0] * vs[1] * vs[2]
            n = round(nv)
            if abs(nv - n) > 1e-3 * max(1.0, abs(nv)) or n == 0:
                continue
            # N(33 s^2) = 33^3 N(s)^2, so N(V)/33^3 must be a rational square
            if n % (33 ** 3):
                continue
            if not issq(n // 33 ** 3):
                continue
            hits.append((c2, c1, c0, n))
            print(f"  *** x^3 + {c2}x^2 + {c1}x + {c0}: N(V) = {n} = "
                  f"33^3 * {n // 33 ** 3} = 33^3 * {math.isqrt(n // 33**3)}^2  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            if len(hits) >= 12:
                break
        if len(hits) >= 12:
            break
    if len(hits) >= 12:
        break
print(f"{tried} totally real cubics with V totally positive; {len(hits)} pass "
      f"the norm filter for V = 33 . square  [{time.time()-t0:.0f}s]",
      flush=True)
