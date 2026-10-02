"""levelt.py p k j [K]: the p-adic plane at level k, G = Cay((Z/p^k)^2, T), modulo the translations by
p^j * Z * (1, 0) (j >= 1, so no translate of a point is adjacent to it: unit vectors are nonzero mod p).
Vertices (x mod p^j, y mod p^k). Writes a K-colouring CNF (edge pinned, or a triangle when 3 is a square)."""
import sys, time
import numpy as np
from collections import defaultdict
p, k, j = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); K = int(sys.argv[4]) if len(sys.argv) > 4 else 4
m = p ** k; mx = p ** j; t0 = time.time()
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
n = mx * m
X = np.repeat(np.arange(mx, dtype=np.int64), m); Y = np.tile(np.arange(m, dtype=np.int64), mx)
E = []
for i in range(0, len(T), 32):
    tb = T[i:i + 32]
    w = (((X[:, None] + tb[None, :, 0]) % mx) * m + (Y[:, None] + tb[None, :, 1]) % m).ravel()
    a = np.repeat(np.arange(n, dtype=np.int64), len(tb))
    assert np.all(a != w)
    E.append(np.unique(np.minimum(a, w) * n + np.maximum(a, w)))
E = np.unique(np.concatenate(E))
print(f"p={p} k={k}, translations by p^{j}(1,0): {n} vertices, {len(E)} edges  [{time.time() - t0:.1f}s]", flush=True)
s3 = roots.get(3 % m, [])
h2 = pow(2, -1, m)
code = lambda x, y: (x % mx) * m + (y % m)
pin = [code(0, 0), code(1, 0)] + ([code(h2, s3[0] * h2)] if s3 else [])
assert len(set(pin)) == len(pin)
tag = f"t{p}k{k}j{j}K{K}"
with open(f"{tag}.cnf", "w") as f:
    f.write(f"p cnf {K * n} {n + K * len(E) + len(pin)}\n")
    for v in range(n):
        f.write(" ".join(str(K * v + c + 1) for c in range(K)) + " 0\n")
    ea, eb = E // n, E % n
    for c in range(K):
        f.write("".join(f"-{K * a + c + 1} -{K * b + c + 1} 0\n" for a, b in zip(ea.tolist(), eb.tolist())))
    for c, v in enumerate(pin):
        f.write(f"{K * v + c + 1} 0\n")
print(f"wrote {tag}.cnf  [{time.time() - t0:.1f}s]", flush=True)
