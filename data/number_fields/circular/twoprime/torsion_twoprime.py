"""For the torsion characters found by torsion_exact.py (one-prime probe, kappa > 2/7), recompute kappa over the
two-prime group <i, rho, sigma> mod M (sigma = (3+2i)/(3-2i)); exact integer arithmetic.  Requires gcd(M, 65) = 1.
Also: over all orders M <= MMAX prime to 65, list every character with kappa > 2/7 over <i, rho, sigma> other than
c and q (expected: none, by Proposition 7)."""
import sys, numpy as np
MMAX = int(sys.argv[1])
def gmul(x, y, M): return ((x[0]*y[0]-x[1]*y[1]) % M, (x[0]*y[1]+x[1]*y[0]) % M)
def ginv(a, M):
    n = (a[0]*a[0] + a[1]*a[1]) % M; ni = pow(n, -1, M); return ((a[0]*ni) % M, (-a[1]*ni) % M)
found = []; examples = {}
for M in range(2, MMAX + 1):
    if M % 5 == 0 or M % 13 == 0: continue
    rho = gmul((2, 1), ginv((2, -1), M), M); sig = gmul((3, 2), ginv((3, -2), M), M)
    grp = {(1 % M, 0)}; frontier = [(1 % M, 0)]
    while frontier:
        nf = []
        for x in frontier:
            for g in (rho, sig, (0, 1 % M)):
                y = gmul(x, g, M)
                if y not in grp: grp.add(y); nf.append(y)
        frontier = nf
    A, B = np.meshgrid(np.arange(M, dtype=np.int64), np.arange(M, dtype=np.int64), indexing="ij")
    kap = np.full((M, M), M, dtype=np.int64)
    for (p, q) in grp:
        v = (A * p + B * q) % M
        np.minimum(kap, np.minimum(v, M - v), out=kap)
    for a, b in zip(*np.nonzero(7 * kap > 2 * M)):
        vals = {int((a*p + b*q) % M) for (p, q) in grp}
        if M % 2 == 0 and vals == {M // 2}: continue
        if M % 3 == 0 and vals <= {M // 3, 2 * M // 3}: continue
        found.append((M, int(a), int(b)))
    if M in (41, 76, 82, 123):
        examples[M] = (len(grp), int(kap[12, 12]) if M == 41 else None)
print("two-prime group sizes / kappa*M of (12+12i)/41:", examples)
print(f"orders M <= {MMAX} prime to 65 with a character of kappa > 2/7 over <i, rho, sigma> other than c, q:", found if found else "none")
