"""Independent Monte Carlo (exact arithmetic) check of Proposition D1: random points of S_1^(r) in [0,N)^2; those that
lie in S_N^(r) (brute-force G_N) must lie in the claimed union (tested by the defining inequalities, not polygons)."""
import random, sys
from fractions import Fraction as Fr
from hgeom import funcs, ev, fl
from sn import GN_brute, in_S, cstar, QN

def in_union(c, k, r):
    N = 5 ** k; s = Fr(1, 2) - r
    # c-component: x = c - c*_N - N m small, all |f(x)| <= s
    T = cstar(N)
    x = (c[0] - T[0] - N * round((c[0] - T[0]) / N), c[1] - T[1] - N * round((c[1] - T[1]) / N))
    if all(abs(ev(f, x)) <= s for j in range(-k, k + 1) for f in funcs(j)):
        return "c"
    for q in QN(N):
        y = (c[0] - q[0] - N * round((c[0] - q[0]) / N), c[1] - q[1] - N * round((c[1] - q[1]) / N))
        ok = True
        for j in range(-k, k + 1):
            for f in funcs(j):
                val = ev(f, q); n = fl(val)
                t = val + ev(f, y)
                if not (n + r <= t <= n + 1 - r):
                    ok = False; break
            if not ok: break
        if ok:
            return "q"
    return None

random.seed(12345)
for (k, r, M) in ((1, Fr(759, 2500), 60000), (2, Fr(759, 2500), 150000), (1, Fr(31, 100), 40000), (2, Fr(1, 3) - Fr(1, 1000), 100000)):
    N = 5 ** k; G = GN_brute(N); D = 10 ** 6
    hits = {"c": 0, "q": 0, None: 0}
    for _ in range(M):
        a, b = random.randrange(N), random.randrange(N)
        t1, t2 = Fr(random.randrange(D + 1), D), Fr(random.randrange(D + 1), D)
        c = (a + r + t1 * (1 - 2 * r), b + r + t2 * (1 - 2 * r))
        if in_S(c, G, r):
            hits[in_union(c, k, r)] += 1
    print(f"k={k} r={r}: {M} samples in S_1; in S_N: c-type {hits['c']}, q-type {hits['q']}, OUTSIDE claimed union {hits[None]}")
    sys.stdout.flush()
