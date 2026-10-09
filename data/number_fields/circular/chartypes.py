"""Map kappa(c) = min over gamma in G_N of ||Re(conj(c) gamma)|| over the torus C / N Z[i] (N = 5^k), on a grid, and
list the local maxima with kappa above a threshold: these are the 'types' of characters that keep every rational
rotation with denominator dividing N inside [r, 1 - r]."""
import sys, numpy as np
from fractions import Fraction as Fr
k = int(sys.argv[1]); h = float(sys.argv[2]); thr = float(sys.argv[3])
N = 5 ** k
rho = complex(3, 4) / 5
G = []
for j in range(-k, k + 1):
    g = rho ** j
    for e in (1, 1j, -1, -1j):
        G.append(e * g)
G = np.array(G)
xs = np.arange(0, N, h)
best = []
for x0 in range(0, len(xs), 200):
    xc = xs[x0:x0 + 200]
    X, Y = np.meshgrid(xc, xs, indexing="ij")
    C = X + 1j * Y
    kap = np.full(C.shape, 0.5)
    for g in G:
        v = (np.conj(C) * g).real
        d = np.abs(v - np.round(v))
        np.minimum(kap, d, out=kap)
    idx = np.argwhere(kap >= thr)
    for a, b in idx:
        best.append((float(kap[a, b]), float(X[a, b]), float(Y[a, b])))
best.sort(reverse=True)
print(f"N = {N}, grid {h}, points with kappa >= {thr}: {len(best)}")
# cluster greedily (points within 3h of an existing cluster centre join it)
clusters = []
for kv, x, y in best:
    for c in clusters:
        dx = (x - c[1] + N / 2) % N - N / 2; dy = (y - c[2] + N / 2) % N - N / 2
        if dx * dx + dy * dy < (0.6) ** 2:
            c[3] += 1; break
    else:
        clusters.append([kv, x, y, 1])
print(f"{len(clusters)} clusters; top by kappa (centre c/N mod 1, as fractions with small denominator):")
for kv, x, y, cnt in sorted(clusters, reverse=True)[:40]:
    fx, fy = Fr(x / N).limit_denominator(60), Fr(y / N).limit_denominator(60)
    print(f"  kappa {kv:.4f}  c = ({x:.3f}, {y:.3f})  c/N ~ ({fx}, {fy})  points {cnt}")
