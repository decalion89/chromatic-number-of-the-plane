"""At which denominator does the realizable direction set start blocking?

The two demands can only be met at once by steps whose denominators stay
bounded: unbounded denominators mean the ball never folds back on itself and
the graph is a tree, while too few steps mean the directions cannot cover
PG(r-1,5).  So the decisive number is the smallest denominator bound d at
which the modulus-one elements of Q(zeta_7) with denominator dividing d
already block every coset 5-colouring.  Below that bound no graph over these
steps can be blocked however large it grows; at or above it, a graph that
realizes all those steps is blocked no matter how it folds.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism

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


# Hilbert 90 over a box: u = alpha / conj(alpha) has modulus one, and every
# modulus-one element of the field arises this way.
found = {}
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
    if d % 5 == 0:          # not 5-integral: no reduction mod 5 exists
        continue
    found.setdefault(d, set()).add(u)

print(f"{sum(len(v) for v in found.values())} modulus-one elements, "
      f"{len(found)} distinct denominators\n")


def direction(u, d):
    """The projective direction of u mod 5, as a primitive integer vector."""
    v = tuple(int(x * d) for x in u)
    g = 0
    for t in v:
        g = gcd(g, abs(t))
    w = tuple(t // g for t in v) if g > 1 else v
    return min(w, tuple(-t for t in w))


print(f"{'denom<=':>8} {'elements':>9} {'directions':>11}  blocks Z/5?")
print("-" * 46)
acc, best = set(), None
for d in sorted(found):
    for u in found[d]:
        acc.add(direction(u, denom(u)))
    S = sorted(acc)
    phi, why = has_homomorphism(S, 5)
    flag = "NO" if phi else "YES  <-- blocks"
    print(f"{d:>8} {sum(len(found[e]) for e in found if e <= d):>9} "
          f"{len(S):>11}  {flag}", flush=True)
    if not phi and best is None:
        best = (d, len(S), sorted(acc))
if best:
    d, ns, S = best
    print(f"\nsmallest blocking denominator bound: {d}  ({ns} directions)")
    import json
    json.dump([list(v) for v in S], open("/tmp/claude-0/-home-user-darwin-50/"
              "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/blocking.json", "w"))
else:
    print("\nno bound in this box blocks")
