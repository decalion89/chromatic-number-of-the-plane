"""Referee E: the 5-adic points c_k = N h + rho^k/5 - 5^(k-1)(2-i) of the paper (Remark (2), Section 9) have
min margin 3/10 over G_N (one-prime probe) but every translate c_k + N m (m mod 13) has min margin < 3/10 over
G(k,1).  Exact."""
from fractions import Fraction as Fr
from math import floor
from ref_rot import rotations_from_generators
for k in range(1, 6):
    N = 5 ** k
    D, rots = rotations_from_generators(k, 1)       # numerators over D = 13 N
    rk = (Fr(1), Fr(0))
    for _ in range(k):
        rk = (rk[0] * Fr(3, 5) - rk[1] * Fr(4, 5), rk[0] * Fr(4, 5) + rk[1] * Fr(3, 5))
    c = (Fr(N, 2) + rk[0] / 5 - 2 * 5 ** (k - 1), Fr(N, 2) + rk[1] / 5 + 5 ** (k - 1))
    def margin(C, keys):
        mn = Fr(1)
        for key in keys:
            p, q = rots[key]
            v = (p * C[0] + q * C[1]) / D
            f = v - floor(v)
            mn = min(mn, f, 1 - f)
        return mn
    keys0 = [key for key in rots if key[2] == 0]
    m0 = margin(c, keys0)
    best = max(margin((c[0] + N * a, c[1] + N * b), list(rots)) for a in range(13) for b in range(13))
    print("k=%d: min margin of c_k over G_N = %s ; max over m of min margin of c_k+Nm over G(k,1) = %s" % (k, m0, best))
