"""minimize_k.py -- shrink a non-4-colourable point set with kissat (+ drat-trim cores) and tabucol.

usage: python3 minimize_k.py PTS.npy TAG [kissat_tlimit_per_call] [max_core_rounds] [deletion_budget_s]
Stage 1 (cores): kissat UNSAT with a DRAT proof (in /dev/shm), drat-trim -c gives the clause core; the vertices
  whose at-least-one clause is in the core (plus the fixed triangle) induce a non-4-colourable subgraph; 4-core;
  repeat until no progress.
Stage 2 (deletion): for each vertex v (lowest degree first) not yet known to be necessary: test S - v
  (tabucol warm-started from the last colouring, then kissat with a time limit). UNSAT -> delete v (4-core).
  SAT with colouring c -> v is necessary; moreover for each colour k with exactly one neighbour w of v coloured k,
  c + (v:=k) colours S - w, so w is necessary too (necessity is inherited by subsets).
Every improvement is saved to min_TAG_pts.npy (+ state json with the triangle position).
"""
import sys, os, time, json, subprocess
import numpy as np
from flat import *
from tabu import tabucol

SC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the scratchpad directory
KISSAT, DRAT = f"{SC}/kissat/build/kissat", f"{SC}/drat-trim/drat-trim"
D, CD, U = directions()
P = np.load(sys.argv[1]).astype(np.int64)
tag = sys.argv[2]
TLK = int(sys.argv[3]) if len(sys.argv) > 3 else 120
MAXCORE = int(sys.argv[4]) if len(sys.argv) > 4 else 8
DELBUDGET = float(sys.argv[5]) if len(sys.argv) > 5 else 1e9
logf = open(f"min_{tag}.log", "a")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); logf.write(s + "\n"); logf.flush()


n = len(P)
E, J = build_edges(P, U)
adj = [set() for _ in range(n)]
for a, b in E:
    adj[a].add(int(b)); adj[b].add(int(a))
stf = sys.argv[1].replace("_pts.npy", "_state.json")
tri = None
if os.path.exists(stf):
    st = json.load(open(stf))
    tri = st.get("triangle") or st.get("triangle_in_saved")
if tri is None:
    tri = find_triangle(n, E.tolist())
tri = tuple(int(x) for x in tri)
assert tri[1] in adj[tri[0]] and tri[2] in adj[tri[0]] and tri[2] in adj[tri[1]]
triset = set(tri)


def core4(S):
    S = set(S) | triset
    stack = [v for v in S if v not in triset and len(adj[v] & S) < 4]
    while stack:
        v = stack.pop()
        if v not in S:
            continue
        S.discard(v)
        for w in adj[v]:
            if w in S and w not in triset and len(adj[w] & S) < 4:
                stack.append(w)
    return S


def sub(S):
    L = sorted(S)
    pos = {v: i for i, v in enumerate(L)}
    Es = [(pos[a], pos[b]) for a in L for b in adj[a] if b in pos and a < b]
    return L, pos, np.array(Es, dtype=np.int64).reshape(-1, 2)


def write_cnf(path, S):
    L, pos, Es = sub(S)
    m = len(L)
    cls = cnf_clauses(m, Es.tolist(), [pos[v] for v in tri])
    with open(path, "w") as f:
        f.write(f"p cnf {4 * m} {len(cls)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in cls))
    return L, pos, Es


def kissat(S, tl, proof=None):
    path = f"/dev/shm/min_{tag}.cnf"
    L, pos, Es = write_cnf(path, S)
    args = ["nice", "-n", "19", KISSAT, f"--time={tl}", path] + ([proof] if proof else [])
    t = time.time()
    out = subprocess.run(args, capture_output=True, text=True).stdout
    dt = time.time() - t
    lines = out.splitlines()
    st = next((l for l in lines if l.startswith("s ")), "s UNKNOWN")
    if "UNSATISFIABLE" in st:
        return False, None, dt, (path, L)
    if "SATISFIABLE" in st:
        lits = set(int(x) for l in lines if l.startswith("v ") for x in l.split()[1:] if int(x) > 0)
        col = {v: next(c for c in range(4) if 1 + 4 * i + c in lits) for i, v in enumerate(L)}
        os.unlink(path)
        return True, col, dt, None
    os.unlink(path)
    return None, None, dt, None


def drat_core(path, L, proof):
    corep = f"/dev/shm/min_{tag}.core"
    t = time.time()
    out = subprocess.run(["nice", "-n", "19", DRAT, path, proof, "-c", corep, "-t", "3000"],
                         capture_output=True, text=True).stdout
    dt = time.time() - t
    ok = "s VERIFIED" in out
    S = set()
    if ok:
        for line in open(corep):
            if line.startswith(("p", "c")):
                continue
            lits = [int(x) for x in line.split()[:-1]]
            if len(lits) == 4 and all(l > 0 for l in lits):
                i = (lits[0] - 1) // 4
                S.add(L[i])
    for f in (corep, proof, path):
        if os.path.exists(f):
            os.unlink(f)
    return ok, S, dt


def save(S, stage):
    L = sorted(S)
    np.save(f"min_{tag}_pts.npy", P[L])
    json.dump({"n": len(L), "stage": stage, "triangle": [L.index(v) for v in tri]}, open(f"min_{tag}_state.json", "w"))


T0 = time.time()
S = core4(range(n))
log(f"=== minimize_k {sys.argv[1]} tag={tag}: {n} vertices, {len(E)} edges; 4-core {len(S)}; triangle {tri}")
# stage 1
for rnd in range(1, MAXCORE + 1):
    proof = f"/dev/shm/min_{tag}.drat"
    r, col, dt, info = kissat(S, 100000, proof)
    if r is not False:
        log(f"core round {rnd}: kissat did not return UNSAT ({r}) -- stop"); sys.exit(1)
    psz = os.path.getsize(proof) / 2 ** 20
    ok, C, dt2 = drat_core(info[0], info[1], proof)
    if not ok:
        log(f"core round {rnd}: drat-trim failed to verify ({dt2:.0f}s)"); break
    C = core4(C)
    log(f"core round {rnd}: kissat UNSAT {dt:.1f}s (proof {psz:.0f} MiB), drat-trim {dt2:.1f}s: core {len(C)} of {len(S)}  [total {time.time() - T0:.0f}s]")
    if len(C) >= len(S):
        break
    S = C
    save(S, f"core{rnd}")
save(S, "cores")
# stage 2
necessary = set(triset)
lastcol = None


def mark(col, v):
    """col colours S - v properly: v is necessary; also every w that is the unique neighbour of v with some colour"""
    necessary.add(v)
    cnt = {}
    for w in adj[v]:
        if w in S:
            if w not in col:       # a neighbour peeled by the 4-core: skip the unique-neighbour rule (safe)
                return
            cnt.setdefault(col[w], []).append(w)
    for k in range(4):
        ws = cnt.get(k, [])
        if len(ws) == 1:
            necessary.add(ws[0])


T2 = time.time()
passes = 0
while True:
    passes += 1
    changed = False
    order = sorted((v for v in S if v not in necessary), key=lambda v: (len(adj[v] & S), v))
    for v in order:
        if v not in S or v in necessary:
            continue
        if time.time() - T2 > DELBUDGET:
            break
        cand = core4(S - {v})
        L, pos, Es = sub(cand)
        init = -np.ones(len(L), dtype=np.int64)
        if lastcol is not None:
            for i, w in enumerate(L):
                init[i] = lastcol.get(w, -1)
        t = time.time()
        ct = tabucol(len(L), Es, init=init, fixed=[(pos[w], i) for i, w in enumerate(tri)], maxiter=1_500_000, seed=v)
        if ct is not None:
            col = {w: int(ct[i]) for i, w in enumerate(L)}
            lastcol = col
            mark(col, v)
            continue
        r, col, dt, info = kissat(cand, TLK)
        if info:
            os.unlink(info[0])
        if r is True:
            lastcol = col
            mark(col, v)
        elif r is False:
            S = cand
            changed = True
            save(S, "deletion")
            log(f"  delete {v}: UNSAT ({dt:.1f}s) -> {len(S)}  [necessary known {len(necessary & S)}; total {time.time() - T0:.0f}s]")
        else:
            necessary.add(v)       # time-out: keep (not proven necessary)
            log(f"  delete {v}: time-out ({dt:.0f}s), kept")
    log(f"deletion pass {passes}: {len(S)} vertices, necessary {len(necessary & S)}  [total {time.time() - T0:.0f}s]")
    if not changed or time.time() - T2 > DELBUDGET:
        break
save(S, "final")
log(f"final: {len(S)} vertices")
