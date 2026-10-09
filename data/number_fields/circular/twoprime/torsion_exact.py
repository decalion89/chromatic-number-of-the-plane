"""(B) Torsion characters of A = Z[1/5][i] of order M prime to 5, for the one-prime probe G_inf = {i^a rho^j}:
c in (1/M)Z[i]/Z[i]; rho acts on Z[i]/M as multiplication by g = (2+i)(2-i)^{-1}; the values at rho^j are periodic.
kappa(c) = min over a period of the distances to Z of Re and Im of conj(c) g^j.  Exact integer arithmetic: the
criterion kappa > 2/7 is 7 * d > 2 * M for the integer distance d = M * kappa.  Lists, for every M <= MMAX, all
characters with kappa > 2/7 other than the bipartite one (all values 1/2) and the four of type q (values 1/3, 2/3).
Also lists, separately, the best kappa among the others (to see how close to 2/7 they come)."""
import sys, numpy as np
MMAX = int(sys.argv[1])
def gmul(x, y, M):
    return ((x[0] * y[0] - x[1] * y[1]) % M, (x[0] * y[1] + x[1] * y[0]) % M)
found = []
best_overall = []
for M in range(2, MMAX + 1):
    if M % 5 == 0: continue
    n = 5 % M; ninv = pow(n, -1, M)             # (2-i)^{-1} = (2+i)/5
    inv = ((2 * ninv) % M, (1 * ninv) % M)
    g = gmul((2, 1), inv, M)
    powers = [(1 % M, 0)]
    while True:
        nxt = gmul(powers[-1], g, M)
        if nxt == (1 % M, 0): break
        powers.append(nxt)
    A, B = np.meshgrid(np.arange(M, dtype=np.int64), np.arange(M, dtype=np.int64), indexing="ij")
    kap = np.full((M, M), M, dtype=np.int64)
    for (p, q) in powers:
        for v in ((A * p + B * q) % M, (A * q - B * p) % M):
            np.minimum(kap, np.minimum(v, M - v), out=kap)
    # exclude the main types exactly: bipartite (M even, all values M/2) and type q (3 | M, values M/3, 2M/3)
    for a, b in zip(*np.nonzero(7 * kap > 2 * M)):
        vals = set()
        for (p, q) in powers:
            vals.add((a * p + b * q) % M); vals.add((a * q - b * p) % M)
        if M % 2 == 0 and vals == {M // 2}: continue
        if M % 3 == 0 and vals <= {M // 3, 2 * M // 3}: continue
        found.append((M, int(a), int(b), int(kap[a, b])))
    # best among non-main (for information)
    mask = np.ones((M, M), dtype=bool)
    if M % 2 == 0: mask[M // 2, M // 2] = False
    if M % 3 == 0:
        for a in (M // 3, 2 * M // 3):
            for b in (M // 3, 2 * M // 3):
                mask[a, b] = False
    mask[0, 0] = False
    kk = np.where(mask, kap, -1)
    i = np.unravel_index(np.argmax(kk), kk.shape)
    best_overall.append((kk[i] / M, M, int(kk[i]), len(powers)))
print(f"orders M <= {MMAX} (prime to 5): characters with kappa > 2/7 other than c and q: {found if found else 'none'}")
best_overall.sort(reverse=True)
print("largest kappa among the other characters:", [(f'{k}/{M}', f'period {P}') for v, M, k, P in best_overall[:8]])
