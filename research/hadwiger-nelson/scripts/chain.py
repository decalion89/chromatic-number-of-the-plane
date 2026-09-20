"""A 4-chromatic graph over a field that can block, with no spindle in it.

The theorem says the critical core must live in a field of residue degree at
least 3 over its real subfield.  The Moser spindle's field is degree 4, so a
spindle can never be that core -- and every 4-critical unit-distance graph
anyone has is multiquadratic, hence out.  Something new is needed.

Chain rhombi instead of pairing them.  A rhombus 0, w, zeta_6 w,
w(1 + zeta_6) forces its tip to the colour of its apex at three colours, and
|1 + zeta_6| = sqrt3 always.  So chaining rhombi in directions w_1, .., w_k
forces every partial sum (1 + zeta_6)(w_1 + .. + w_j) to the colour of the
origin, and if

    | w_1 + .. + w_k | = 1/sqrt3

the last one is at distance 1 from the origin -- forced equal and adjacent, so
no three-colouring exists.  k = 2 is the Moser spindle and needs
u.vbar + ubar.v = -5/3, i.e. (-5 + sqrt-11)/6, which Q(zeta_21) provably does
not contain.  k = 3 asks instead for

    r(u,v) + r(u,w) + r(v,w) = -4/3,     r(a,b) = a.bbar + abar.b,

three terms of the real subfield summing to a rational -- a real search, not a
closed form, and one nothing rules out.
"""
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField

K = CycloField(21)
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


def den(t):
    d = 1
    for x in t:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


Z6 = K.neg(K.mul(K.zeta(7), K.zeta(7)))
steps = set()
for coeffs in itertools.product(range(-2, 3), repeat=6):
    if not any(coeffs):
        continue
    a, p = K.zero(), ONE
    for c in coeffs:
        a = K.add(a, tuple(Fraction(c) * x for x in p))
        p = K.mul(p, K.zeta(1) if False else K.zeta(21))
    try:
        u = K.mul(a, inv(K.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if K.norm2(u) != ONE:
        continue
    v = u
    for _ in range(6):
        steps.add(v)
        v = K.mul(v, Z6)
steps = sorted(steps)
dens = {}
for u in steps:
    dens[den(u)] = dens.get(den(u), 0) + 1
print(f"{len(steps)} unit steps of Q(zeta_21); denominators "
      + ", ".join(f"{k}x{v}" for k, v in sorted(dens.items())[:12]), flush=True)

t0 = time.time()
TARGET2 = tuple([Fraction(-5, 3)] + [Fraction(0)] * (D - 1))
TARGET3 = tuple([Fraction(-4, 3)] + [Fraction(0)] * (D - 1))
r = {}
for i, u in enumerate(steps):
    for j in range(i + 1, len(steps)):
        v = steps[j]
        t = K.mul(u, K.conj(v))
        r[(i, j)] = K.add(t, K.conj(t))
print(f"{len(r)} pair values precomputed  [{time.time()-t0:.0f}s]", flush=True)

pairs = [k for k, v in r.items() if v == TARGET2]
print(f"k = 2 (the Moser spindle): {len(pairs)} solutions "
      "-- the no-spindle theorem says this must be 0", flush=True)

found, n = [], len(steps)
for i in range(n):
    for j in range(i + 1, n):
        rij = r[(i, j)]
        for k in range(j + 1, n):
            s = K.add(K.add(rij, r[(i, k)]), r[(j, k)])
            if s == TARGET3:
                found.append((i, j, k))
                if len(found) <= 3:
                    print(f"  *** k = 3 CHAIN FOUND: steps {i}, {j}, {k}  "
                          f"[{time.time()-t0:.0f}s]", flush=True)
    if i % 20 == 0:
        print(f"  ... i = {i}/{n}, {len(found)} so far "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"k = 3: {len(found)} chains  [{time.time()-t0:.0f}s]", flush=True)
