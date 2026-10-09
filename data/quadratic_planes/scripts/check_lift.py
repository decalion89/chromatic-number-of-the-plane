"""check_lift.py p k IN.npy: an independent check that the preimage at level k of the level-(k-1) points in IN.npy
has no proper 4-colouring in G_k = Cay((Z/p^k)^2, T_k). Written separately from liftcore.py:
- T_k from the rational parametrisation of the circle and Hensel lifting is NOT used; instead every (a, b) mod p^k
  is tested directly, column by column;
- the vertex set is listed point by point (x0 + p^(k-1) i, y0 + p^(k-1) j) and indexed by a dict;
- the colour-major encoding (variable c * n + v + 1) and a pinned triangle found by a fresh search;
then kissat with a DRAT proof and drat-trim."""
import sys, os, subprocess
import numpy as np
p, k, IN = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # kissat/ and drat-trim/ live here
m, mp = p ** k, p ** (k - 1)
T = []
a = np.arange(m, dtype=np.int64)
for b in range(m):
    ok = np.nonzero((a * a + b * b) % m == 1)[0]
    T.extend((int(x), b) for x in ok)
T = sorted(set(T))
assert len(T) == (p + 1) * p ** (k - 1), len(T)
if IN.endswith(".json"):                 # padic11.json: the points of level k - 1 as [x, y]
    import json
    prev = [x * mp + y for x, y in json.load(open(IN))[f"level{k - 1}"]]
else:
    prev = [int(c) for c in np.load(IN)]
pts = sorted({((c // mp + mp * i) % m, (c % mp + mp * j) % m) for c in prev for i in range(p) for j in range(p)})
idx = {q: n for n, q in enumerate(pts)}
n = len(pts)
E = set()
for v, (x, y) in enumerate(pts):
    for (s, t) in T:
        w = idx.get(((x + s) % m, (y + t) % m))
        if w is not None and w != v:
            E.add((min(v, w), max(v, w)))
E = sorted(E)
print(f"check: p={p} level {k}: {n} points, {len(E)} edges, {len(T)} unit vectors", flush=True)
adj = {}
for u, v in E:
    adj.setdefault(u, set()).add(v); adj.setdefault(v, set()).add(u)
tri = None
for u in range(n - 1, -1, -1):          # search from the other end this time
    for v in adj.get(u, ()):
        c = adj[u] & adj[v]
        if c:
            tri = (u, v, min(c)); break
    if tri:
        break
var = lambda v, c: c * n + v + 1
path = f"/dev/shm/checklift_{os.getpid()}"
with open(path + ".cnf", "w") as f:
    pin = tri if tri else E[0]
    f.write(f"p cnf {4 * n} {n + 4 * len(E) + len(pin)}\n")
    for v in range(n):
        f.write(" ".join(str(var(v, c)) for c in range(4)) + " 0\n")
    for c in range(4):
        f.write("".join(f"-{var(u, c)} -{var(v, c)} 0\n" for u, v in E))
    for c, v in enumerate(pin):
        f.write(f"{var(v, c)} 0\n")
print(f"pinned {'triangle' if tri else 'edge'} {pin}", flush=True)
import hashlib
digest = hashlib.sha256(open(path + ".cnf", "rb").read()).hexdigest()
print(f"formula sha256 {digest}", flush=True)
k_out = subprocess.run([f"{SC}/kissat/build/kissat", path + ".cnf", path + ".drat"], capture_output=True, text=True).stdout
print("kissat:", [l for l in k_out.splitlines() if l.startswith("s ")], flush=True)
d_out = subprocess.run([f"{SC}/drat-trim/drat-trim", path + ".cnf", path + ".drat", "-t", "20000"], capture_output=True, text=True).stdout
print("drat-trim:", [l for l in d_out.splitlines() if l.startswith("s ")], flush=True)
LOG = os.environ.get("LOGDIR")
if LOG:
    os.makedirs(LOG, exist_ok=True)
    open(f"{LOG}/level{k}.kissat.log", "w").write(k_out)
    open(f"{LOG}/level{k}.drat-trim.log", "w").write(d_out)
os.unlink(path + ".cnf")
os.unlink(path + ".drat")
