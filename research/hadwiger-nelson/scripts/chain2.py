"""Chained rhombi over Q(zeta_21): a 4-chromatic graph with no spindle in it.

A rhombus 0, w, zeta_6 w, w(1 + zeta_6) forces its tip to the apex's colour at
three colours, and |1 + zeta_6| = sqrt3 whatever w is.  Chain them in
directions w_1, .., w_k and every partial sum (1 + zeta_6)(w_1 + .. + w_j) is
forced to the origin's colour; if |w_1 + .. + w_k| = 1/sqrt3 the last one is
also ADJACENT to the origin, and no three-colouring exists.

k = 2 is the Moser spindle and needs (-5 + sqrt-11)/6, which Q(zeta_21) does
not contain -- checked here as a control.  k = 3 asks for

    |u + v + w|^2 = 1/3,   i.e.   r(u,v) + r(u,w) + r(v,w) = -4/3

with r(a,b) = a.bbar + abar.b in the real subfield: three irrational terms
summing to a rational.  Nothing rules it out, and a graph found this way lives
over Q(zeta_21), which the residue-degree theorem says CAN block, and contains
no spindle, which is what stops the critical core from collapsing.
"""
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField
from hn.homcol import denominator_29_directions

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


Z6 = K.neg(K.mul(K.zeta(7), K.zeta(7)))
steps = set()
for coeffs in itertools.product(range(-2, 3), repeat=6):
    if not any(coeffs):
        continue
    a, p = K.zero(), ONE
    for c in coeffs:
        a = K.add(a, tuple(Fraction(c) * x for x in p))
        p = K.mul(p, K.zeta(1))
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
    d = 1
    for x in u:
        d = d * x.denominator // gcd(d, x.denominator)
    dens[d] = dens.get(d, 0) + 1
print(f"{len(steps)} unit steps of Q(zeta_21); denominators "
      + ", ".join(f"{k}x{v}" for k, v in sorted(dens.items())[:14]), flush=True)

t0 = time.time()
SC = 3 * 29 * 29
n = len(steps)
r = [[None] * n for _ in range(n)]
bad = 0
for i in range(n):
    ci = K.conj(steps[i])
    for j in range(i + 1, n):
        t = K.mul(steps[j], ci)
        val = K.add(t, K.conj(t))
        try:
            iv = tuple(int(x * SC) if (x * SC).denominator == 1 else None
                       for x in val)
        except Exception:
            iv = None
        if iv is None or any(z is None for z in iv):
            bad += 1
            iv = None
        r[i][j] = r[j][i] = iv
print(f"{n*(n-1)//2} pair values, {bad} off the common denominator "
      f"[{time.time()-t0:.0f}s]", flush=True)

T2 = tuple([-5 * SC // 3] + [0] * (D - 1))
T3 = tuple([-4 * SC // 3] + [0] * (D - 1))
two = [(i, j) for i in range(n) for j in range(i + 1, n) if r[i][j] == T2]
print(f"k = 2 control (the Moser spindle): {len(two)} solutions "
      "-- the no-spindle theorem says 0", flush=True)

found = []
for i in range(n):
    ri = r[i]
    for j in range(i + 1, n):
        rij = ri[j]
        if rij is None:
            continue
        rj = r[j]
        need = tuple(a - b for a, b in zip(T3, rij))
        for k in range(j + 1, n):
            a, b = ri[k], rj[k]
            if a is None or b is None:
                continue
            if all(x + y == z for x, y, z in zip(a, b, need)):
                found.append((i, j, k))
                print(f"  *** k = 3 CHAIN: steps {i}, {j}, {k}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    if i % 25 == 0:
        print(f"  ... i = {i}/{n}, {len(found)} so far "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"k = 3: {len(found)} chains total  [{time.time()-t0:.0f}s]", flush=True)
