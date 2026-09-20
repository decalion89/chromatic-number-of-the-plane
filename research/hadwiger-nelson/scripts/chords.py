"""Which chords do the BLOCKING rotations offer?

A spindle's rotation rho must satisfy |1 - rho| = 1/d, where d is a distance
carrying a forced-same pair.  If rho is itself a blocking step, its direction
becomes a genuine edge of the construction rather than a pendant -- which is
what "load bearing" means here.

So the question is what values |1 - rho|^2 takes over the denominator-29
modulus-one elements of Q(zeta_7).  Each one names a distance d = 1/|1 - rho|
that a spindle could use, and only the RATIONAL values are usable, since a
chord in a degree-3 real subfield will not match a lattice distance.

Note one candidate in advance. rho with rational real part must generate a
quadratic field, and the quadratic subfield of Q(zeta_7) is Q(sqrt-7), so
rho = (a + b sqrt-7)/c with a^2 + 7b^2 = c^2. At c = 29 that has the solution
a = 27, b = 4, giving |1 - rho|^2 = 2 - 54/29 = 4/29 and d = sqrt(29)/2.
"""
import sys, itertools
from fractions import Fraction
from math import gcd, isqrt
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField

F = CycloField(7)
D = F.degree
one = F.rational(1)


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def denom(u):
    d = 1
    for x in u:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


steps = {}
for coeffs in itertools.product(range(-3, 4), repeat=4):
    a = tuple(Fraction(c) for c in coeffs) + (Fraction(0), Fraction(0))
    if not any(a):
        continue
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) != one:
        continue
    d = denom(u)
    if d % 5 == 0 or d > 29:
        continue
    chord = F.norm2(F.sub(one, u))
    steps.setdefault(d, []).append((u, chord))

print("chords |1 - rho|^2 by denominator, and which are rational:")
for d in sorted(steps):
    rats = sorted({c[0] for _, c in steps[d] if not any(c[1:])})
    print(f"  denominator {d:3}: {len(steps[d])} steps, rational chords "
          f"{[str(r) for r in rats] if rats else 'none'}", flush=True)

print("\nusable spindle distances d = 1/|1 - rho| from the rational chords:")
seen = set()
for d in sorted(steps):
    for u, c in steps[d]:
        if any(c[1:]):
            continue
        q = c[0]
        if q == 0 or q in seen:
            continue
        seen.add(q)
        dsq = Fraction(1) / q
        print(f"  denominator {d}: chord^2 = {q}, so d^2 = {dsq}"
              f"   (d = {float(dsq) ** 0.5:.5f})", flush=True)
