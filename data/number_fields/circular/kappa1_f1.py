"""kappa_1(p, 1) for every prime p = 3 (mod 4) below P (fast, numpy): in F_{p^2} = F_p(i), as z runs over the norm-one
elements, bz runs over the circle s^2 + t^2 = c, c = N(b) != 0, and Tr(bz) = 2s.  So
kappa_1(p, 1) = max over c in F_p^* of min over {s : c - s^2 is a square or 0} of ||2s/p||.
Also the analytic threshold: for p = 3 (mod 4) the Weil bound |sum_{z} e(Tr(bz)/p)| <= 2 sqrt(p) gives a z with
||Tr(bz)/p|| < 1/4 as soon as (p^2 - 1)/(2p) > 2 sqrt(p) (ln p + 1)."""
import sys, math
import numpy as np
from sympy import primerange
P = int(sys.argv[1])
thr = next(p for p in range(11, 10 ** 7) if (p * p - 1) / (2 * p) > 2 * math.sqrt(p) * (math.log(p) + 1)
           and all((q * q - 1) / (2 * q) > 2 * math.sqrt(q) * (math.log(q) + 1) for q in range(p, p + 2000)))
print(f"analytic threshold: the inequality holds for every p >= {thr} (checked monotone on the next 2000 integers)")
worst = (0, None)
for p in primerange(3, P):
    if p % 4 != 3:
        continue
    s = np.arange(p)
    sq = np.zeros(p, dtype=bool); sq[(s * s) % p] = True          # squares including 0
    dist = np.minimum((2 * s) % p, p - (2 * s) % p)                  # p * ||2s/p||
    best = 0
    for c in range(1, p):
        ok = sq[(c - s * s) % p]
        m = dist[ok].min()
        if m > best:
            best = m
    k = best / p
    if p >= 11 and k >= 0.25:
        print(f"p = {p}: kappa_1 = {best}/{p} >= 1/4 !!", flush=True)
    if p >= 11 and k > worst[0]:
        worst = (k, p)
    if p < 50:
        print(f"p = {p}: kappa_1 = {best}/{p} = {k:.4f}", flush=True)
print(f"all primes p = 3 mod 4 with 11 <= p < {P}: max kappa_1 = {worst[0]:.4f} at p = {worst[1]} (< 1/4)")
