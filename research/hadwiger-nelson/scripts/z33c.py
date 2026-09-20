"""Q(zeta_33) done properly: build sqrt-11 instead of hoping to stumble on it.

The first pass enumerated a/abar over a box of low powers of zeta_33 and found
no Moser spindle at all -- which is impossible, since sqrt-11 lies in
Q(zeta_11) and rho = (5 + sqrt-11)/6 has modulus one.  The box was the fault:
the Gauss sum

    g = sum_{k=1..10} (k|11) zeta_11^k,    g^2 = -11

needs zeta_33^{3k} up to k = 10, far outside a box of the first seven powers.
Recorded as a caught error rather than quietly widened: an enumeration that
misses an element it was built to find will report a clean zero and look like
a theorem.

Built directly here, and the unit steps are generated MULTIPLICATIVELY -- they
form a group, so products of a few generators reach far further than any box.
"""
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism

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


Z11 = K.zeta(3)
Z3 = K.zeta(11)
Z6 = K.neg(K.mul(Z3, Z3))
QR = {1, 3, 4, 5, 9}
g = K.zero()
p = ONE
for k in range(1, 11):
    p = K.mul(p, Z11) if k > 1 else Z11
    g = K.add(g, p if k in QR else K.neg(p))
print(f"Gauss sum g: g^2 = -11? {K.mul(g, g) == K.rational(-11)}", flush=True)
RHO = tuple(Fraction(1, 6) * x for x in K.add(K.rational(5), g))
print(f"rho = (5 + sqrt-11)/6: |rho|^2 = 1? {K.norm2(RHO) == ONE}; "
      f"rho + rhobar = 5/3? "
      f"{K.add(RHO, K.conj(RHO)) == K.rational(Fraction(5, 3))}", flush=True)

t0 = time.time()
extra = set()
for coeffs in itertools.product(range(-1, 2), repeat=5):
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
print(f"{len(extra)} unit steps from Q(zeta_11)  [{time.time()-t0:.0f}s]",
      flush=True)

zp, z = [], ONE
for _ in range(6):
    zp.append(z)
    z = K.mul(z, Z6)
gens = zp + [RHO, K.conj(RHO), K.neg(RHO), K.neg(K.conj(RHO))] \
    + sorted(extra)[:40]
steps, frontier = {ONE}, [ONE]
for _ in range(2):
    nxt = []
    for x in frontier:
        for gq in gens:
            y = K.mul(x, gq)
            if y not in steps:
                steps.add(y)
                nxt.append(y)
    frontier = nxt
    if len(steps) > 6000:
        break
steps = sorted(steps)
print(f"{len(steps)} unit steps after closing under products  "
      f"[{time.time()-t0:.0f}s]", flush=True)

TWO = K.rational(Fraction(-5, 3))
FOUR3 = K.rational(Fraction(-4, 3))
re = {i: K.add(u, K.conj(u)) for i, u in enumerate(steps)}
sp = [i for i in range(len(steps)) if re[i] == TWO]
print(f"k = 2 (Moser spindles): {len(sp)} rotations -- must be nonzero now",
      flush=True)

n = len(steps)
found = []
for i in range(n):
    ci = K.conj(steps[i])
    ri = re[i]
    for j in range(i + 1, n):
        t = K.mul(steps[j], ci)
        if K.add(K.add(ri, re[j]), K.add(t, K.conj(t))) == FOUR3:
            found.append((i, j))
            if len(found) <= 3:
                print(f"  *** k = 3 CHAIN: a = {i}, b = {j}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    if i % 100 == 0:
        print(f"  ... {i}/{n}, {len(found)} chains [{time.time()-t0:.0f}s]",
              flush=True)
print(f"k = 3: {len(found)} chains  [{time.time()-t0:.0f}s]", flush=True)
