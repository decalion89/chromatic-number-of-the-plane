"""liftcore.py p k OUT [IN.npy] [K]: obstruction lifting in the p-adic plane (p = 3 mod 4).
The vertex set is all of (Z/p^k)^2, or, with IN.npy (codes x * p^(k-1) + y of points at level k - 1), the
preimage of that set at level k. The graph is the induced subgraph of G_k = Cay((Z/p^k)^2, T_k). Asks kissat
whether it is K-colourable (K = 4; a triangle pinned if 3 is a square mod p, else an edge); if not, shrinks to the
drat-trim core (vertices whose 'some colour' clause the proof uses), repeating while it shrinks, and saves the
codes of the core (level k) to OUT.npy. A preimage that is not K-colourable proves that G_k is not K-colourable."""
import sys, os, subprocess, time
import numpy as np
from collections import defaultdict
p, k, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
IN = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] != "-" else None
K = int(sys.argv[5]) if len(sys.argv) > 5 else 4
SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # kissat/ and drat-trim/ live here
KISSAT, DRAT = f"{SC}/kissat/build/kissat", f"{SC}/drat-trim/drat-trim"
m = p ** k; t0 = time.time()
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
SAME = os.environ.get("SAME") == "1"     # IN holds level-k codes: use them as they are
if IN and SAME:
    V = np.unique(np.load(IN).astype(np.int64))
elif IN:
    mp = p ** (k - 1)
    if IN.endswith(".json"):             # padic11.json: the points of level k - 1 as [x, y]
        import json
        prev = np.array([x * mp + y for x, y in json.load(open(IN))[f"level{k - 1}"]], dtype=np.int64)
    else:
        prev = np.load(IN).astype(np.int64)
    px, py = prev // mp, prev % mp
    lifts = np.arange(p, dtype=np.int64) * mp
    X = (px[:, None, None] + lifts[None, :, None] + 0 * lifts[None, None, :]).ravel()
    Y = (py[:, None, None] + 0 * lifts[None, :, None] + lifts[None, None, :]).ravel()
    V = np.unique(X * m + Y)
else:
    V = np.arange(m * m, dtype=np.int64)
vx, vy = V // m, V % m


def lookup(codes):
    """index in V of each code, or -1 (V is sorted)"""
    i = np.searchsorted(V, codes)
    i[i >= len(V)] = 0
    return np.where(V[i] == codes, i, -1)


E = []
CH = max(1, 4_000_000 // len(V))
for i in range(0, len(T), CH):
    tb = T[i:i + CH]
    w = lookup((((vx[:, None] + tb[None, :, 0]) % m) * m + (vy[:, None] + tb[None, :, 1]) % m).ravel())
    a = np.repeat(np.arange(len(V), dtype=np.int64), len(tb))
    ok = w >= 0
    a, w = a[ok], w[ok]
    E.append(np.unique(np.minimum(a, w) * len(V) + np.maximum(a, w)))
E = np.unique(np.concatenate(E)); n = len(V)
EA, EB = E // n, E % n
print(f"p={p} level {k}{' (preimage of ' + IN + ')' if IN else ''}: {n} vertices, {len(E)} edges  [{time.time() - t0:.1f}s]", flush=True)
s3 = roots.get(3 % m, [])
TMP = f"/dev/shm/liftcore_{os.getpid()}"


def solve(S):
    """S: array of vertex indices. Returns (sat?, core vertex indices or None)."""
    Sset = np.zeros(n, dtype=bool); Sset[S] = True
    keep = Sset[EA] & Sset[EB]
    ea, eb = EA[keep], EB[keep]
    loc = -np.ones(n, dtype=np.int64); loc[S] = np.arange(len(S))
    la, lb = loc[ea], loc[eb]
    # pin: a triangle if possible, else an edge
    adj = defaultdict(set)
    for a_, b_ in zip(la[:200000].tolist(), lb[:200000].tolist()):
        adj[a_].add(b_); adj[b_].add(a_)
    pin = None
    if s3:
        for a_ in list(adj)[:2000]:
            for b_ in adj[a_]:
                c_ = adj[a_] & adj[b_]
                if c_:
                    pin = [a_, b_, next(iter(c_))]; break
            if pin: break
    if pin is None:
        pin = [int(la[0]), int(lb[0])]
    with open(TMP + ".cnf", "w") as f:
        f.write(f"p cnf {K * len(S)} {len(S) + K * len(la) + len(pin)}\n")
        for v in range(len(S)):
            f.write(" ".join(str(K * v + c + 1) for c in range(K)) + " 0\n")
        for c in range(K):
            f.write("".join(f"-{K * x + c + 1} -{K * y + c + 1} 0\n" for x, y in zip(la.tolist(), lb.tolist())))
        for c, v in enumerate(pin):
            f.write(f"{K * v + c + 1} 0\n")
    r = subprocess.run([KISSAT, "--time=3000", TMP + ".cnf", TMP + ".drat"], capture_output=True, text=True).returncode
    if r == 10:
        return True, None
    if r != 20:
        return None, None
    o = subprocess.run([DRAT, TMP + ".cnf", TMP + ".drat", "-c", TMP + ".core", "-t", "5000"], capture_output=True, text=True).stdout
    os.unlink(TMP + ".drat")
    if "s VERIFIED" not in o:
        return None, None
    C = set(pin)
    for line in open(TMP + ".core"):
        if line[0] in "pc":
            continue
        t = [int(x) for x in line.split()[:-1]]
        if len(t) == K and all(x > 0 for x in t):
            C.add((t[0] - 1) // K)
    return False, S[np.array(sorted(C), dtype=np.int64)]


S = np.arange(n, dtype=np.int64)
r, C = solve(S)
print(f"  all {n}: {'K-colourable' if r else ('NOT ' + str(K) + '-colourable (DRAT verified)' if r is False else 'unknown')}  [{time.time() - t0:.0f}s]", flush=True)
while r is False and len(C) < len(S):
    S = C
    np.save(out + ".npy", V[S])            # every core found is itself not K-colourable: keep it
    r, C = solve(S)
    print(f"  core {len(S)}: {'NOT colourable, next core ' + str(len(C)) if r is False else r}  [{time.time() - t0:.0f}s]", flush=True)
if r is False:
    np.save(out + ".npy", V[S])
    print(f"saved {out}.npy: {len(S)} level-{k} points with no proper {K}-colouring", flush=True)
for ext in (".cnf", ".core"):
    if os.path.exists(TMP + ext):
        os.unlink(TMP + ext)
