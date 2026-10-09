# Sanity check of the remark on Cohen's conjecture: in Q(sqrt(-d)) inside C, every unit
# z = (a + e sqrt(-d))/D (a^2 + d e^2 = D^2, gcd = 1) has p not dividing D for an odd prime p | d,
# and a/D = +-1 mod p; so colouring x + y sqrt(-d) (x, y p-integral) by a proper colouring of the
# p-cycle at (x mod p) is proper on every edge.
from math import gcd, isqrt
for d in (3, 7, 11, 15, 19, 23, 35, 47, 59, 71, 83, 131, 6, 10, 14):
    p = next(q for q in range(3, d + 1, 2) if d % q == 0)
    cnt = 0
    for D in range(1, 2000):
        for e in range(0, isqrt(D * D // d) + 1):
            a2 = D * D - d * e * e
            a = isqrt(a2)
            if a * a != a2 or gcd(gcd(a, e), D) != 1:
                continue
            cnt += 1
            assert D % p != 0, (d, a, e, D)
            r = a * pow(D, -1, p) % p
            assert r in (1, p - 1), (d, a, e, D, r)
    print(f"d={d}: p={p}, {cnt} primitive units (a,e >= 0, D < 2000), all = +-1 mod p")
