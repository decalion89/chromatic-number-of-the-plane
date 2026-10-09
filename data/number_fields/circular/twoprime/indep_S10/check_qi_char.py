"""Line 1020-1021: the 7-torsion points (65/7)(a+bi) as characters c -> Re(conj(c) gamma) of Q(i) at rotations whose
denominator does not divide 65; and the correct 7-adic character gamma -> lambda(gamma mod 7)/7."""
from fractions import Fraction as Fr
import math

c = (Fr(65, 7), Fr(325, 7))           # (65/7)(1+5i), the point of lambda(a+bi) = 2a+3b
best = None
for n in range(2, 200):
    for a in range(-n, n + 1):
        b2 = n * n - a * a
        b = math.isqrt(b2)
        if b * b != b2 or math.gcd(math.gcd(a, b), n) != 1:
            continue
        for bb in (b, -b):
            g = (Fr(a, n), Fr(bb, n))
            w = (c[0] * g[0] + c[1] * g[1], c[0] * g[1] - c[1] * g[0])   # conj(c) g
            for x in w:
                f = x - math.floor(x)
                m = min(f, 1 - f)
                if best is None or m < best[0]:
                    best = (m, n, a, bb)
    if best and best[0] < Fr(1, 20):
        break
print("smallest margin of Re/Im(conj((65/7)(1+5i)) gamma) over rotations found:", best,
      "(float %.4f); 2/7 = %.4f" % (float(best[0]), 2 / 7))
# the 7-adic character itself at that rotation
m, n, a, b = best
inv = pow(n, -1, 7)
r = ((a * inv) % 7, (b * inv) % 7)
for (p, q) in (r, (r[1], (-r[0]) % 7)):        # gamma and -i gamma
    print("  lambda(gamma~)/7 at gamma, -i gamma:", Fr((2 * p + 3 * q) % 7, 7))
