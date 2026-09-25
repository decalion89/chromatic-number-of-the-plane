"""Search the curve directly in r, not through the enumerated steps.

The closing condition depends only on r = a + abar:

    y1^2 = 3(4 - r^2)      a exists as a unit step with that real part
    y2^2 = 3(8 - 12r - 9r^2)   its closing partner b lies in K

Both are conditions on ONE element of the real subfield F, so the natural
search is over F itself rather than over pairs of enumerated steps -- which is
what the first sweeps did, and why every chain they returned had r of degree
2, sitting inside Q(sqrt33) and dragging the whole graph back into the
spindle's dead degree-4 field.

r is built from the basis eta_j = zeta^j + zeta^-j of F = Q(zeta_33)+, degree
10, with small rational coefficients.  Squareness is tested by reduction at
primes 1 mod 33, where Q(zeta_33) splits completely so reduction is evaluation
at a root of unity; a single non-residue rejects r.  The embeddings are
checked too: |sigma(r)| <= 2 at every real place, since sigma(a)sigma(abar) = 1
holds at all of them in a CM field.

What is wanted is a survivor with deg Q(r) >= 3.  That is what lifts the
chain's module out of rank 4.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, time, cmath
from fractions import Fraction
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.cyclotomic import CycloField

K = CycloField(33)
D = K.degree
ONE = K.rational(1)
Z = K.zeta(1)

# eta_j = zeta^j + zeta^-j, j coprime to 33, one per conjugate pair
reps, seen = [], set()
for j in range(1, 33):
    if __import__("math").gcd(j, 33) != 1:
        continue
    if (33 - j) in seen:
        continue
    seen.add(j)
    reps.append(j)
ETA = []
for j in reps:
    zj = ONE
    for _ in range(j):
        zj = K.mul(zj, Z)
    ETA.append(K.add(zj, K.conj(zj)))
print(f"F = Q(zeta_33)+ : {len(ETA)} basis elements eta_j, j in {reps}",
      flush=True)

PRIMES = [67, 199, 331, 397, 463, 661, 727]
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


EMB = []
for k in range(1, 33):
    if __import__("math").gcd(k, 33) == 1:
        EMB.append(k)


def embeddings(v):
    return [sum(complex(float(x)) * cmath.exp(2j * cmath.pi * i * k / 33)
                for i, x in enumerate(v)) for k in EMB]


def deg_of(v):
    """Degree of Q(v): the number of distinct embeddings of v."""
    vals = embeddings(v)
    out = []
    for z in vals:
        if not any(abs(z - w) < 1e-7 for w in out):
            out.append(z)
    return len(out)


t0 = time.time()
THREE = K.rational(3)
found = []
DENOMS = (1, 2, 3, 6)
RANGE = range(-3, 4)
count = 0
for dn in DENOMS:
    for coeffs in itertools.product(RANGE, repeat=4):
        if not any(coeffs):
            continue
        r = K.zero()
        for c, e in zip(coeffs, ETA):
            if c:
                r = K.add(r, tuple(Fraction(c, dn) * x for x in e))
        count += 1
        w1 = K.mul(THREE, K.sub(K.rational(4), K.mul(r, r)))
        if not square_mod(w1):
            continue
        w2 = K.mul(THREE, K.sub(K.sub(K.rational(8), tuple(12 * x for x in r)),
                                tuple(9 * x for x in K.mul(r, r))))
        if not square_mod(w2):
            continue
        if max(abs(z.real) for z in embeddings(r)) > 2 + 1e-9:
            continue
        d = deg_of(r)
        found.append((d, coeffs, dn))
        print(f"  *** r = {coeffs}/{dn}: BOTH squares, deg Q(r) = {d}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"  denominator {dn} done, {count} tried, {len(found)} found "
          f"[{time.time()-t0:.0f}s]", flush=True)
hi = [f for f in found if f[0] >= 3]
print(f"{len(found)} points on the curve, {len(hi)} with deg Q(r) >= 3  "
      f"[{time.time()-t0:.0f}s]", flush=True)
