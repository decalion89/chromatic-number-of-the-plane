"""levelp.py p k SUB [K]: the p-adic unit-distance plane at level k (p = 3 mod 4), G = Cay((Z/p^k)^2, T),
T = {x^2 + y^2 = 1 mod p^k}, modulo the rotations H = {t in T : t = 1 mod p^SUB} (SUB = k: no symmetry).
No H-orbit contains an edge, since (h - 1) z = 0 mod p is never a unit vector. Writes a K-colouring CNF
(K = 4 by default) of the orbit graph, with a triangle pinned if 3 is a square mod p and an edge otherwise."""
import sys, time
import numpy as np
from collections import defaultdict
p, k, SUB = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); K = int(sys.argv[4]) if len(sys.argv) > 4 else 4
assert p % 4 == 3
m = p ** k; t0 = time.time()
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
assert len(T) == (p + 1) * p ** (k - 1), len(T)
q = p ** SUB
H = T[(T[:, 0] % q == 1 % q) & (T[:, 1] % q == 0)]
print(f"p={p} k={k}: |T| = {len(T)}, |H| = {len(H)}", flush=True)
X = np.repeat(np.arange(m, dtype=np.int64), m); Y = np.tile(np.arange(m, dtype=np.int64), m)
rep = X * m + Y
for a, b in H.tolist():
    np.minimum(rep, ((X * a - Y * b) % m) * m + (X * b + Y * a) % m, out=rep)
del X, Y
reps, inv = np.unique(rep, return_inverse=True)
n = len(reps)
rx, ry = reps // m, reps % m
E = []
for j in range(0, len(T), 32):
    tb = T[j:j + 32]
    w = inv[(((rx[:, None] + tb[None, :, 0]) % m) * m + (ry[:, None] + tb[None, :, 1]) % m).ravel()]
    a = np.repeat(np.arange(n, dtype=np.int64), len(tb))
    assert np.all(a != w), "an orbit contains an edge"
    E.append(np.unique(np.minimum(a, w) * n + np.maximum(a, w)))
E = np.unique(np.concatenate(E))
print(f"orbit graph: {n} vertices, {len(E)} edges  [{time.time() - t0:.1f}s]", flush=True)
s3 = roots.get(3 % m, [])
h2 = pow(2, -1, m)
if s3:
    pin = [inv[0], inv[1 * m + 0], inv[(h2 % m) * m + (s3[0] * h2) % m]]
else:
    pin = [inv[0], inv[1 * m + 0]]
assert len(set(pin)) == len(pin)
tag = f"p{p}k{k}s{SUB}K{K}"
with open(f"{tag}.cnf", "w") as f:
    f.write(f"p cnf {K * n} {n + K * len(E) + len(pin)}\n")
    for v in range(n):
        f.write(" ".join(str(K * v + c + 1) for c in range(K)) + " 0\n")
    ea, eb = E // n, E % n
    for c in range(K):
        f.write("".join(f"-{K * x + c + 1} -{K * y + c + 1} 0\n" for x, y in zip(ea.tolist(), eb.tolist())))
    for c, v in enumerate(pin):
        f.write(f"{K * v + c + 1} 0\n")
np.save(f"{tag}_reps.npy", reps); np.save(f"{tag}_edges.npy", E)
print(f"wrote {tag}.cnf ({'triangle' if s3 else 'edge'} pinned)  [{time.time() - t0:.1f}s]", flush=True)
