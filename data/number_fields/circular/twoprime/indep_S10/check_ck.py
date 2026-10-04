"""Remark (2) of Section 9: c_k = N h + rho^k/5 - 5^(k-1)(2-i) lies in S_N^{3/10} with least margin exactly 3/10, and
its distances to N E_c and N E_q are at least sqrt5 5^(k-1) - 1/5 and (7/6) 5^(k-1) - 1/5."""
from fractions import Fraction as Fr
import math


def cdiv(x, y):
    n = y[0] * y[0] + y[1] * y[1]
    return (Fr(x[0] * y[0] + x[1] * y[1], n), Fr(x[1] * y[0] - x[0] * y[1], n))


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z, e = cdiv((1, 0), z), -e
    for _ in range(e):
        w = cmul(w, z)
    return w


RHO = cdiv((2, 1), (2, -1))
for k in range(1, 9):
    N = 5 ** k
    rk = cpow(RHO, k)
    c = (Fr(N, 2) + rk[0] / 5 - 2 * 5 ** (k - 1), Fr(N, 2) + rk[1] / 5 + 5 ** (k - 1))
    m = Fr(1)
    for j in range(-k, k + 1):
        w = cmul((c[0], -c[1]), cpow(RHO, j))
        for x in w:
            f = x - math.floor(x)
            m = min(m, f, 1 - f)
    # distance^2 to N E_c = N h + N Z[i] and to N E_q = N(a+bi)/3 + N Z[i]: minimise over nearby lattice points
    def d2(centre):
        best = None
        for u in range(-3, 4):
            for v in range(-3, 4):
                p = (centre[0] + N * u, centre[1] + N * v)
                dd = (c[0] - p[0]) ** 2 + (c[1] - p[1]) ** 2
                best = dd if best is None or dd < best else best
        return best
    dc = d2((Fr(N, 2), Fr(N, 2)))
    dq = min(d2((Fr(N * a, 3), Fr(N * b, 3))) for a in (1, 2) for b in (1, 2))
    lc = math.sqrt(5) * 5 ** (k - 1) - 0.2
    lq = 7 / 6 * 5 ** (k - 1) - 0.2
    print(f"k={k}: least margin over G_N = {m}; dist to NE_c = {math.sqrt(dc):.4f} >= {lc:.4f}: {math.sqrt(dc) >= lc}; "
          f"dist to NE_q = {math.sqrt(dq):.4f} >= {lq:.4f}: {math.sqrt(dq) >= lq}")
