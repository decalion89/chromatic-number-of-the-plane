"""Referee check 8: field facts.
 - sqrt3 in Q(sqrt d) iff 3d is a square (d squarefree, d != 3): not for d = 7, 31, 15, 39.
 - a unit triangle in F^2 forces sqrt3 in F: with u = B - A, v = C - A, |u|^2 = |v|^2 = |u - v|^2 = 1 gives
   u.v = 1/2 and det(u, v)^2 = |u|^2 |v|^2 - (u.v)^2 = 3/4 (identity checked symbolically below on random data).
 - Theorem B / Corollary B6 conditions for real quadratic fields as stated in the paper: (a) iff d != 3 mod 4,
   (b) iff d != 2 mod 3.
 - sqrt7 and sqrt31 in Q_3: Hensel lifting of x^2 = d modulo 3^k.
"""
import math, random
from fractions import Fraction as Fr

for d in (7, 31, 15, 39):
    s = math.isqrt(3 * d)
    print(f'd = {d}: 3d = {3 * d} square: {s * s == 3 * d} -> sqrt3 in Q(sqrt{d}): {s * s == 3 * d};'
          f' (a) [d != 3 mod 4]: {d % 4 != 3}; (b) [d != 2 mod 3]: {d % 3 != 2};'
          f' so chi = {2 if d % 4 != 3 else (3 if d % 3 != 2 else ">=4")}')

# Lagrange identity det^2 = |u|^2|v|^2 - (u.v)^2 on random rationals
ok = True
for _ in range(1000):
    u = (Fr(random.randint(-50, 50), random.randint(1, 50)), Fr(random.randint(-50, 50), random.randint(1, 50)))
    v = (Fr(random.randint(-50, 50), random.randint(1, 50)), Fr(random.randint(-50, 50), random.randint(1, 50)))
    det = u[0] * v[1] - u[1] * v[0]; dot = u[0] * v[0] + u[1] * v[1]
    ok &= det * det == (u[0] ** 2 + u[1] ** 2) * (v[0] ** 2 + v[1] ** 2) - dot * dot
print('Lagrange identity det^2 = |u|^2|v|^2 - (u.v)^2 holds on 1000 random pairs:', ok,
      '-> unit triangle gives det^2 = 1 - 1/4 = 3/4, so sqrt3 = 2 det in F')


def hensel_sqrt(a, p, k):
    x = next(t for t in range(1, p) if (t * t - a) % p == 0)
    m = p
    for _ in range(k - 1):
        m2 = m * p
        # Newton step x <- x - (x^2 - a)/(2x) mod m2
        x = (x - (x * x - a) * pow(2 * x, -1, m2)) % m2
        m = m2
    return x, m


for a in (7, 31):
    x, m = hensel_sqrt(a, 3, 40)
    print(f'sqrt{a} in Z_3: x = {x} mod 3^40, x^2 - {a} divisible by 3^40: {(x * x - a) % m == 0}')
print('3 is not a square in Q_3 (odd valuation), so Q_3 has no sqrt3 and its plane no unit triangle')
