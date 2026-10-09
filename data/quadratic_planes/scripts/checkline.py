"""checkline.py p k IN.npy|all [samples]: checks the line description used by liftline.py against the definition.
For every pair of level-k blocks z, w of the set at unit distance mod p^k (or a random sample of them), and every
pair of next digits x, y: z~ + p^k x and w~ + p^k y are at unit distance mod p^(k+1)  <=>  ubar.(y - x) = sign*c,
with ubar, sign, c computed as in liftline.py.  Also checks that no two points of one block are adjacent, and
that blocks at non-unit distance mod p^k have no edges between them (sampled)."""
import sys, random
import numpy as np
from collections import defaultdict
p, k, IN = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; samples = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
m = p ** k; M = m * p; inv2 = pow(2, -1, p)
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = [(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])]
S = list(range(m * m)) if IN == "all" else sorted(set(np.load(IN).tolist()))
Sset = set(S)
pairs = []
for z in S:
    for (a, b) in T:
        w = ((z // m + a) % m) * m + (z % m + b) % m
        if w in Sset and w > z:
            pairs.append((z, w))
random.seed(1)
if len(pairs) > samples:
    pairs = random.sample(pairs, samples)
adjM = lambda P, Q: ((Q[0] - P[0]) ** 2 + (Q[1] - P[1]) ** 2 - 1) % M == 0
bad = 0
for z, w in pairs:
    zx, zy, wx, wy = z // m, z % m, w // m, w % m
    dx, dy = wx - zx, wy - zy
    e = ((dx * dx + dy * dy - 1) // m) % p; c = (-e * inv2) % p
    cu, nu = (dx % p) * p + dy % p, ((-dx) % p) * p + (-dy) % p
    ub = min(cu, nu); sg = 1 if cu <= nu else -1; ux, uy = ub // p, ub % p
    for x0 in range(p):
        for x1 in range(p):
            for y0 in range(p):
                for y1 in range(p):
                    P = (zx + m * x0, zy + m * x1); Q = (wx + m * y0, wy + m * y1)
                    model = (ux * (y0 - x0) + uy * (y1 - x1)) % p == (sg * c) % p
                    if adjM(P, Q) != model:
                        bad += 1
# no edges inside a block, none between blocks at non-unit distance (sampled)
for _ in range(20000):
    z, w = random.choice(S), random.choice(S)
    x, y = (random.randrange(p), random.randrange(p)), (random.randrange(p), random.randrange(p))
    P = (z // m + m * x[0], z % m + m * x[1]); Q = (w // m + m * y[0], w % m + m * y[1])
    unit_k = ((w // m - z // m) ** 2 + (w % m - z % m) ** 2 - 1) % m == 0 and z != w
    if adjM(P, Q) and not unit_k:
        bad += 1
print(f"p={p} k={k}: checked {len(pairs)} block edges x {p ** 4} digit pairs, + 20000 random pairs: {bad} mismatches")
