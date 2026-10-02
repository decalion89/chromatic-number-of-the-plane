"""quotlevel.py p k vx vy [K]: K-colourings of the level-k p-adic plane G_k invariant under the translations by
p^(k-1) * (vx, vy) (a subgroup of order p; no point is adjacent to its translates, since unit vectors are
nonzero mod p). Builds the quotient graph and asks kissat."""
import sys, os, subprocess, time
import numpy as np
from collections import defaultdict
p, k, vx, vy = (int(a) for a in sys.argv[1:5]); K = int(sys.argv[5]) if len(sys.argv) > 5 else 4
SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # kissat/ and drat-trim/ live here
m = p ** k; n = m * m; t0 = time.time(); h = p ** (k - 1)
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
X = np.arange(n, dtype=np.int64) // m; Y = np.arange(n, dtype=np.int64) % m
canon = np.full(n, n, dtype=np.int64)
for j in range(p):
    canon = np.minimum(canon, ((X + j * h * vx) % m) * m + (Y + j * h * vy) % m)
reps, cls = np.unique(canon, return_inverse=True); nc = len(reps)
rx, ry = reps // m, reps % m
E = set()
for (a, b) in T:
    w = cls[(((rx + a) % m) * m + (ry + b) % m)]
    assert np.all(w != np.arange(nc))
    lo = np.minimum(np.arange(nc), w); hi = np.maximum(np.arange(nc), w)
    E.update((lo * nc + hi).tolist())
E = np.array(sorted(E), dtype=np.int64); ea, eb = E // nc, E % nc
print(f"p={p} k={k} mod <{h}*({vx},{vy})>: {nc} classes, {len(E)} edges  [{time.time() - t0:.1f}s]", flush=True)
cnf = f"quot_{p}_{k}_{vx}_{vy}_K{K}.cnf"
with open(cnf, "w") as f:
    f.write(f"p cnf {K * nc} {nc + K * len(E) + 2}\n")
    for v in range(nc):
        f.write(" ".join(str(K * v + c + 1) for c in range(K)) + " 0\n")
    for c in range(K):
        f.write("".join(f"-{K * x + c + 1} -{K * y + c + 1} 0\n" for x, y in zip(ea.tolist(), eb.tolist())))
    f.write(f"{K * ea[0] + 1} 0\n{K * eb[0] + 2} 0\n")
r = subprocess.run([f"{SC}/kissat/build/kissat", "--time=3000", cnf, cnf + ".drat"], capture_output=True, text=True)
res = [l for l in r.stdout.splitlines() if l.startswith("s ")]
print("  kissat:", res, f"[{time.time() - t0:.0f}s]", flush=True)
if res == ["s UNSATISFIABLE"]:
    o = subprocess.run([f"{SC}/drat-trim/drat-trim", cnf, cnf + ".drat", "-t", "20000"], capture_output=True, text=True).stdout
    print("  drat-trim:", [l for l in o.splitlines() if l.startswith("s ")], f"[{time.time() - t0:.0f}s]", flush=True)
    open(cnf + ".kissat.log", "w").write(r.stdout); open(cnf + ".drat-trim.log", "w").write(o)
if os.path.exists(cnf + ".drat"):
    os.unlink(cnf + ".drat")
