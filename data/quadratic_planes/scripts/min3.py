"""min3.py -- shrink a non-3-colourable unit-distance graph over Q(sqrt47) to a vertex-critical one.

usage: python3 min3.py STATE.json PTS.npy OUT_PREFIX
1. drat-trim cores: kissat UNSAT with a DRAT proof, the vertices whose at-least-one clause is in the clause core,
   the 3-core of the induced subgraph (a vertex of degree <= 2 never decides 3-colourability); repeat while it shrinks.
2. deletion: for each vertex, lowest degree first, delete it if the rest is still not 3-colourable (kissat);
   otherwise it is necessary (kept). The result is vertex-critical: every vertex is necessary.
Writes OUT_PREFIX_pts.npy and OUT_PREFIX.json (points, D, U, edges)."""
import sys, os, json, subprocess
import numpy as np

SC = os.environ.get("SC", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KISSAT, DRAT = f"{SC}/kissat/build/kissat", f"{SC}/drat-trim/drat-trim"
st = json.load(open(sys.argv[1]))
P = np.load(sys.argv[2]).astype(np.int64)
out = sys.argv[3]
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


def sat(S, proof=None):
    path = f"/dev/shm/min3_{os.getpid()}.cnf"
    L, E = cnf(S, path)
    args = [KISSAT, "--time=600", path] + ([proof] if proof else [])
    o = subprocess.run(args, capture_output=True, text=True).stdout
    s = next((l for l in o.splitlines() if l.startswith("s ")), "s UNKNOWN")
    return ("UNSAT" in s and False) or ("SATISFIABLE" in s and "UNSAT" not in s and True) or (None if "UNSAT" not in s else False), L


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


S = core3(range(n))
print(f"start {n}, 3-core {len(S)}", flush=True)
while True:
    proof = f"/dev/shm/min3_{os.getpid()}.drat"
    r, L = sat(S, proof)
    assert r is False, "not UNSAT?"
    o = subprocess.run([DRAT, f"/dev/shm/min3_{os.getpid()}.cnf", proof, "-c", f"/dev/shm/min3_{os.getpid()}.core", "-t", "3000"],
                       capture_output=True, text=True).stdout
    assert "s VERIFIED" in o
    C = set()
    for line in open(f"/dev/shm/min3_{os.getpid()}.core"):
        if line[0] in "pc":
            continue
        t = [int(x) for x in line.split()[:-1]]
        if len(t) == 3 and all(x > 0 for x in t):
            C.add(L[(t[0] - 1) // 3])
    C = core3(C)
    r2, _ = sat(C)
    if r2 is not False:            # the fixed edge may lie outside: keep its endpoints
        C = core3(C | set(L[:0]))
        r2, _ = sat(C)
    print(f"core round: {len(S)} -> {len(C)} (UNSAT: {r2 is False})", flush=True)
    if r2 is not False or len(C) >= len(S):
        break
    S = C
necessary = set()
changed = True
while changed:
    changed = False
    for v in sorted(S, key=lambda v: (len(adj[v] & S), v)):
        if v in necessary or v not in S:
            continue
        T = core3(S - {v})
        r, _ = sat(T)
        if r is False:
            S = T; changed = True
            print(f"delete {v}: -> {len(S)}", flush=True)
        else:
            necessary.add(v)
print(f"vertex-critical: {len(S)} vertices", flush=True)
L = sorted(S)
pos = {v: k for k, v in enumerate(L)}
E = [[pos[a], pos[b]] for a in L for b in adj[a] if b in pos and a < b]
np.save(out + "_pts.npy", np.array([pts[v] for v in L], dtype=np.int64))
json.dump({"D": D, "U": [list(u) for u in U], "points": [list(pts[v]) for v in L], "edges": E}, open(out + ".json", "w"))
print(f"saved {out}: {len(L)} vertices, {len(E)} edges", flush=True)
