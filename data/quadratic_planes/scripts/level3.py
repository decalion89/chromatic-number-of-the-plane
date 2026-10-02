"""level3.py k [L [S]]: the 11-adic unit-distance plane at level k, G_k = Cay((Z/11^k)^2, T_k), T_k = {x^2 + y^2 = 1},
modulo the rotations by H = {t in T_k : t = 1 mod 11^S} (default S = 1; a colouring constant on H-orbits; no H-orbit
contains an edge, since (h - 1) z is never a unit). With L = 1 also modulo the translations by 11^(k-1) (Z/11) (1, 0)
(no loops either). Writes the orbit graph and a 4-colouring CNF (a unit triangle pinned) to g{k}[t][s{S}].cnf in the
current directory. 2 October: python3 level3.py 3 (15 961 orbits) and python3 level3.py 3 1 2 (27 951 orbits) are
both UNSAT with kissat; python3 level3.py 2 is UNSAT as it must be (G_2 has no 4-colouring)."""
import sys, time
import numpy as np
k = int(sys.argv[1]); TRANS = len(sys.argv) > 2 and sys.argv[2] == "1"
m = 11 ** k
t0 = time.time()
x = np.arange(m, dtype=np.int64)
# T_k: for each x, the y with y^2 = 1 - x^2 mod m
sq = (x * x) % m
from collections import defaultdict
roots = defaultdict(list)
for y, s in enumerate(sq.tolist()):
    roots[s].append(y)
T = [(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])]
T = np.array(T, dtype=np.int64)
print(f"k={k}: |T_k| = {len(T)} (expected {12 * 11 ** (k - 1)})  [{time.time() - t0:.1f}s]", flush=True)
SUB = int(sys.argv[3]) if len(sys.argv) > 3 else 1
H = T[(T[:, 0] % 11 ** SUB == 1) & (T[:, 1] % 11 ** SUB == 0)]
print(f"|H| = {len(H)}", flush=True)
# all points, encoded z = x * m + y
X = np.repeat(np.arange(m, dtype=np.int64), m); Y = np.tile(np.arange(m, dtype=np.int64), m)
code = X * m + Y
rep = code.copy()
for a, b in H.tolist():
    xr = (X * a - Y * b) % m; yr = (X * b + Y * a) % m
    c = xr * m + yr
    np.minimum(rep, c, out=rep)
    if TRANS:
        for j in range(1, 11):
            c2 = ((xr + j * 11 ** (k - 1)) % m) * m + yr
            np.minimum(rep, c2, out=rep)
del X, Y
# with translations, the minimum over the group generated needs closure: iterate rep <- rep[rep] until stable
it = 0
while True:
    nr = rep[rep]
    # also re-apply the generators to the representatives
    if np.array_equal(nr, rep):
        break
    rep = nr; it += 1
reps, inv = np.unique(rep, return_inverse=True)
n = len(reps)
print(f"orbits: {n} (closure iterations {it})  [{time.time() - t0:.1f}s]", flush=True)
# edges: for each orbit representative z and t in T: orbit of z + t
rx, ry = reps // m, reps % m
E = []
for j in range(0, len(T), 64):
    tb = T[j:j + 64]
    nx = (rx[:, None] + tb[None, :, 0]) % m; ny = (ry[:, None] + tb[None, :, 1]) % m
    w = inv[(nx * m + ny).ravel()]
    a = np.repeat(np.arange(n, dtype=np.int64), len(tb))
    lo, hi = np.minimum(a, w), np.maximum(a, w)
    assert np.all(lo != hi), "a loop: an orbit contains an edge"
    E.append(np.unique(lo * n + hi))
E = np.unique(np.concatenate(E))
print(f"orbit graph: {n} vertices, {len(E)} edges  [{time.time() - t0:.1f}s]", flush=True)
# a triangle: 0, (1, 0), (1/2, s/2) with s^2 = 3
s = next(y for y in range(m) if (y * y - 3) % m == 0)
h2 = pow(2, -1, m)
tri = [inv[0], inv[1 * m + 0], inv[(h2 % m) * m + (s * h2) % m]]
print("triangle orbits", tri, flush=True)
tag = f"g{k}" + ("t" if TRANS else "") + (f"s{SUB}" if SUB != 1 else "")
with open(f"{tag}.cnf", "w") as f:
    K = 4
    f.write(f"p cnf {K * n} {n + K * len(E) + 3}\n")
    for v in range(n):
        f.write(" ".join(str(K * v + c + 1) for c in range(K)) + " 0\n")
    ea, eb = E // n, E % n
    for c in range(K):
        f.write("".join(f"-{K * a + c + 1} -{K * b + c + 1} 0\n" for a, b in zip(ea.tolist(), eb.tolist())))
    for c, v in enumerate(tri):
        f.write(f"{K * v + c + 1} 0\n")
np.save(f"{tag}_reps.npy", reps); np.save(f"{tag}_edges.npy", E)
print(f"wrote {tag}.cnf  [{time.time() - t0:.1f}s]", flush=True)
