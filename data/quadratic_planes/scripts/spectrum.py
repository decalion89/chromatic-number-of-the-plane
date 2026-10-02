"""spectrum.py p L: normalised eigenvalues mu_l(N) of the p-adic unit-distance plane at levels l = 1..L
(p = 3 mod 4). A character of (Z/p^k)^2 is p^j * eta with eta primitive mod p^l (l = k - j); its normalised
eigenvalue is the average of e(<eta, u>/p^l) over the unit circle T_l mod p^l, which depends only on
N = |eta|^2 mod p^l (rotations act transitively on primitive vectors of a given norm). Prints, per level, the
minimum and maximum of mu_l and the Hoffman ratio bound -mu/(1 - mu) for the minimum over levels <= l."""
import sys
import numpy as np
from collections import defaultdict
p, L = int(sys.argv[1]), int(sys.argv[2])
best = 1.0
for l in range(1, L + 1):
    m = p ** l
    roots = defaultdict(list)
    for y in range(m):
        roots[(y * y) % m].append(y)
    T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
    # one primitive eta per unit norm N: eta = (a, b) with a^2 + b^2 = N
    etas = {}
    for a in range(m):
        for b in range(min(m, 400)):
            N = (a * a + b * b) % m
            if N % p and N not in etas:
                etas[N] = (a, b)
        if len(etas) == (p - 1) * p ** (l - 1):
            break
    assert len(etas) == (p - 1) * p ** (l - 1), len(etas)
    Ns = sorted(etas)
    E = np.array([etas[N] for N in Ns], dtype=np.int64)
    mu = np.zeros(len(Ns))
    for i in range(0, len(Ns), 256):
        ph = (E[i:i + 256, 0:1] * T[None, :, 0] + E[i:i + 256, 1:2] * T[None, :, 1]) % m
        mu[i:i + 256] = np.cos(2 * np.pi * ph / m).mean(axis=1)
    lo, hi = mu.min(), mu.max()
    best = min(best, lo)
    sq = np.array([pow(N, (p - 1) // 2, p) == 1 for N in Ns])
    print(f"p={p} level {l}: |T| = {len(T)}, {len(Ns)} norms; mu in [{lo:.5f}, {hi:.5f}]; "
          f"N square mod p: [{mu[sq].min():.4f}, {mu[sq].max():.4f}], non-square: [{mu[~sq].min():.4f}, {mu[~sq].max():.4f}]; "
          f"Hoffman ratio over levels <= {l}: {-best / (1 - best):.5f} (1/4 = 0.25, 1/5 = 0.2)", flush=True)
