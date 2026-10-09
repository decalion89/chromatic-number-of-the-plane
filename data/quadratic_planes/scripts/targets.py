"""targets.py d D [kmax]: the vectors m of the module spanned by U_D (sums of at most kmax unit vectors) at a
'spindle distance': 4|m|^2 - 1 = t^2 with t in Q(sqrt d). Prints them by number of steps, then |m|^2."""
import sys
from fractions import Fraction as Fr
from math import isqrt
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from units_fast import units_fast
d, D = int(sys.argv[1]), int(sys.argv[2]); kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 3
U = np.array(units_fast(d, D), dtype=np.int64)

def sqrt_int_K(P0, Q0):
    """integers p, q with (p + q sqrt d)^2 = P0 + Q0 sqrt d, or None (Z[sqrt d] is the ring of integers, d = 3 mod 4)"""
    N = P0 * P0 - d * Q0 * Q0
    if N < 0: return None
    n = isqrt(N)
    if n * n != N: return None
    for s in (n, -n):
        if (P0 + s) % 2: continue
        p2 = (P0 + s) // 2
        if p2 < 0: continue
        p = isqrt(p2)
        if p * p != p2: continue
        for pp in ({p, -p} if p else {0}):
            if pp == 0:
                if Q0 == 0 and P0 % d == 0 and P0 >= 0:
                    q2 = P0 // d; q = isqrt(q2)
                    if q * q == q2: return (0, q)
                continue
            if Q0 % (2 * pp) == 0:
                q = Q0 // (2 * pp)
                if pp * pp + d * q * q == P0: return (pp, q)
    return None

S = {(0, 0, 0, 0): 0}
frontier = [(0, 0, 0, 0)]
found = []
for k in range(1, kmax + 1):
    new = []
    for p in frontier:
        for u in U.tolist():
            q = (p[0] + u[0], p[1] + u[1], p[2] + u[2], p[3] + u[3])
            if q not in S:
                S[q] = k; new.append(q)
    frontier = new
    for q in new:
        a, b, c, e = q
        X = a * a + d * b * b + c * c + d * e * e; Y = 2 * (a * b + c * e)      # |m|^2 = (X + Y sqrt d)/D^2
        r = sqrt_int_K(4 * X - D * D, 4 * Y)                                     # 4|m|^2 - 1 = (P0 + Q0 sqrt d)/D^2
        if r is not None:
            found.append((k, Fr(X, D * D), Fr(Y, D * D), q, (Fr(r[0], D), Fr(r[1], D))))
    print(f"k={k}: {len(new)} new vectors; spindle vectors so far {len(found)}", flush=True)
found.sort(key=lambda f: (f[0], float(f[1]) + float(f[2]) * d ** 0.5))
for k, X, Y, q, t in found[:40]:
    print(f"steps {k}: m = {q}  |m|^2 = {X} + {Y} sqrt{d}  t = {t[0]} + {t[1]} sqrt{d}")
