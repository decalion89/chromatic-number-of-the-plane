"""Torsion characters of Z[1/5][i] (order M prime to 5) are c in (1/M)Z[i]/Z[i]; rho acts on Z[i]/M as multiplication
by g = (2+i)(2-i)^{-1} mod M, so the values at rho^j are periodic in j and
    kappa_inf(c) = min over j in a period, and over the coordinates, of ||(Re, Im)(conj(c) g^j)||   (values mod 1).
Lists, for each M, the best kappa_inf among characters that are neither the bipartite one (all values 1/2) nor of
type q (all values in {1/3, 2/3}).  Exact integer arithmetic (values are integers mod M)."""
import sys, numpy as np
from math import gcd
MMAX = int(sys.argv[1]); THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.29
def inv_mod_gauss(a, b, M):
    # inverse of a+bi in Z[i]/M: (a - bi)/(a^2+b^2) mod M
    n = (a * a + b * b) % M
    ninv = pow(n, -1, M)
    return ((a * ninv) % M, (-b * ninv) % M)
def gmul(x, y, M):
    return ((x[0] * y[0] - x[1] * y[1]) % M, (x[0] * y[1] + x[1] * y[0]) % M)
results = []
for M in range(2, MMAX + 1):
    if M % 5 == 0: continue
    inv = inv_mod_gauss(2, -1, M)
    g = gmul((2, 1), inv, M)
    # period of g
    powers = [(1 % M, 0)]
    while True:
        nxt = gmul(powers[-1], g, M)
        if nxt == (1 % M, 0): break
        powers.append(nxt)
    A, B = np.meshgrid(np.arange(M), np.arange(M), indexing="ij")
    alive = np.ones((M, M), dtype=bool)
    kap = np.full((M, M), M // 2, dtype=np.int64)    # kappa * M as integer distance
    for (p, q) in powers:
        # conj(c) * (p + q i) with c = (A + B i)/M: Re = A p + B q, Im = A q - B p  (mod M)
        for v in ((A * p + B * q) % M, (A * q - B * p) % M):
            d = np.minimum(v, M - v)
            np.minimum(kap, d, out=kap)
    best = None
    for a, b in zip(*np.nonzero(kap > THR * M)):
        val = kap[a, b] / M
        # type tests: values set
        vals = set()
        for (p, q) in powers:
            vals.add(((a * p + b * q) % M, (a * q - b * p) % M))
        flat = {x for pr in vals for x in pr}
        if flat == {M // 2} and M % 2 == 0: continue                       # bipartite c*
        if M % 3 == 0 and flat <= {M // 3, 2 * M // 3}: continue            # type q
        if best is None or val > best[0]:
            best = (val, (a, b), len(powers), sorted(flat))
    if best:
        results.append((best[0], M, best[1], best[2]))
        print(f"M={M}: best non-main kappa_inf = {best[0]:.6f} = {kap[best[1]]}/{M} at c = ({best[1][0]} + {best[1][1]} i)/{M}, period {best[2]}")
        sys.stdout.flush()
results.sort(reverse=True)
print("top:", [(round(v, 6), M) for v, M, c, p in results[:10]])
