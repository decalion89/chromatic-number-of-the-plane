"""min3fast.py -- like min3.py, faster on large graphs.

usage: python3 min3fast.py STATE.json PTS.npy OUT_PREFIX [seeds]
1. drat-trim cores as in min3.py, but each round runs kissat with several random seeds and keeps the
   smallest core whose 3-core is still not 3-colourable; repeat while it shrinks.
2. deletion: for each vertex v, lowest degree first, look for a 3-colouring of the 3-core of S - v with tabu
   search first (a short budget); only if tabu fails, ask kissat. v is deleted when kissat says UNSAT.
   The result is vertex-critical (every remaining vertex was shown necessary by an explicit colouring),
   and certify_q.py checks it again from scratch.
Writes OUT_PREFIX_pts.npy and OUT_PREFIX.json (points, D, U, edges), like min3.py."""
import sys, os, json, subprocess
import numpy as np
from tabu import tabucol

SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KISSAT, DRAT = f"{SC}/kissat/build/kissat", f"{SC}/drat-trim/drat-trim"
st = json.load(open(sys.argv[1]))
P = np.load(sys.argv[2]).astype(np.int64)
out = sys.argv[3]
SEEDS = int(sys.argv[4]) if len(sys.argv) > 4 else 4
D, U = st["D"], [tuple(u) for u in st["U"]]
pts = [tuple(int(x) for x in p) for p in P.tolist()]
idx = {p: i for i, p in enumerate(pts)}
n = len(pts)
adj = [set() for _ in range(n)]
for i, p in enumerate(pts):
    for u in U:
        j = idx.get(tuple(a + b for a, b in zip(p, u)))
        if j is not None:
            adj[i].add(j); adj[j].add(i)
TMP = f"/dev/shm/min3f_{os.getpid()}"


def cnf(S, path):
    L = sorted(S)
    pos = {v: k for k, v in enumerate(L)}
    E = [(pos[a], pos[b]) for a in L for b in adj[a] if b in pos and a < b]
    cl = [[3 * k + c + 1 for c in range(3)] for k in range(len(L))]
    cl += [[-(3 * a + c + 1), -(3 * b + c + 1)] for a, b in E for c in range(3)]
    a0, b0 = E[0]
    cl += [[3 * a0 + 1], [3 * b0 + 2]]
    with open(path, "w") as f:
        f.write(f"p cnf {3 * len(L)} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    return L, E


def kissat(S, proof=None, seed=0):
    L, E = cnf(S, TMP + ".cnf")
    args = [KISSAT, "--time=900", f"--seed={seed}", TMP + ".cnf"] + ([proof] if proof else [])
    r = subprocess.run(args, capture_output=True, text=True).returncode
    return {10: True, 20: False}.get(r), L


def core3(S):
    S = set(S)
    stack = [v for v in S if len(adj[v] & S) < 3]
    while stack:
        v = stack.pop()
        if v not in S:
            continue
        S.discard(v)
        for w in adj[v]:
            if w in S and len(adj[w] & S) < 3:
                stack.append(w)
    return S


def drat_core(S, seed):
    proof = TMP + ".drat"
    r, L = kissat(S, proof, seed)
    if r is not False:
        return None
    o = subprocess.run([DRAT, TMP + ".cnf", proof, "-c", TMP + ".core", "-t", "3000"], capture_output=True, text=True).stdout
    os.unlink(proof)
    if "s VERIFIED" not in o:
        return None
    C = set()
    for line in open(TMP + ".core"):
        if line[0] in "pc":
            continue
        t = [int(x) for x in line.split()[:-1]]
        if len(t) == 3 and all(x > 0 for x in t):
            C.add(L[(t[0] - 1) // 3])
    return core3(C)


S = core3(range(n))
print(f"start {n}, 3-core {len(S)}", flush=True)
while True:
    best = None
    for seed in range(SEEDS):
        C = drat_core(S, seed)
        if C is None or len(C) >= len(S):
            continue
        r, _ = kissat(C)
        if r is False and (best is None or len(C) < len(best)):
            best = C
    print(f"core round: {len(S)} -> {len(best) if best else '-'}", flush=True)
    if best is None:
        break
    S = best

necessary = set()
changed = True
tests = tabu_keeps = kissat_calls = 0
while changed:
    changed = False
    for v in sorted(S, key=lambda v: (len(adj[v] & S), v)):
        if v in necessary or v not in S:
            continue
        T = core3(S - {v})
        tests += 1
        L = sorted(T)
        pos = {w: k for k, w in enumerate(L)}
        E = np.array([(pos[a], pos[b]) for a in L for b in adj[a] if b in pos and a < b], dtype=np.int64)
        col = tabucol(len(L), E, k=3, maxiter=300_000, seed=tests) if len(E) else np.zeros(len(L), dtype=np.int64)
        if col is not None:
            necessary.add(v); tabu_keeps += 1
            continue
        kissat_calls += 1
        r, _ = kissat(T)
        if r is False:
            S = T; changed = True
            print(f"delete {v}: -> {len(S)}  [tests {tests}, tabu keeps {tabu_keeps}, kissat {kissat_calls}]", flush=True)
        else:
            necessary.add(v)
print(f"vertex-critical: {len(S)} vertices  [tests {tests}, tabu keeps {tabu_keeps}, kissat {kissat_calls}]", flush=True)
L = sorted(S)
pos = {v: k for k, v in enumerate(L)}
E = [[pos[a], pos[b]] for a in L for b in adj[a] if b in pos and a < b]
np.save(out + "_pts.npy", np.array([pts[v] for v in L], dtype=np.int64))
json.dump({"D": D, "U": [list(u) for u in U], "points": [list(pts[v]) for v in L], "edges": E}, open(out + ".json", "w"))
print(f"saved {out}: {len(L)} vertices, {len(E)} edges", flush=True)
for ext in (".cnf", ".core"):
    if os.path.exists(TMP + ext):
        os.unlink(TMP + ext)
