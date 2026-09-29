"""grow3.py -- colouring-guided growth (tabucol repair first, kissat as the complete solver).

usage: python3 grow3.py SEED TAG [maxn] [kissat_tlimit] [tabu_iters] [K_frac]
  SEED: H, HP, G1, B2, B1U, or PTS.npy (optionally with PTS_state.json holding "triangle" and "col")
Round: colour the graph (tabucol from the previous colouring; if that fails, kissat with a time limit);
candidates y = x + u (x in graph, u in U) not in the graph; blocked = the graph-neighbours of y see all 4 colours;
add the best blocked ones (most neighbours, then in Z[zeta21], then smallest trace norm). Stop at UNSAT
(kissat), maxn, or kissat time-out. State saved every round (grow_TAG_pts.npy, grow_TAG_state.json).
"""
import sys, os, time, json, subprocess
import numpy as np
from flat import *
from seeds import heptagon, g1_points, EMB
from tabu import tabucol
from lns import lns_repair

SC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the scratchpad directory
KISSAT = f"{SC}/kissat/build/kissat"
D, CD, U = directions()
Ua = np.array(U, dtype=np.int64)
seed, tag = sys.argv[1], sys.argv[2]
MAXN = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
TL = int(sys.argv[4]) if len(sys.argv) > 4 else 1200
TABU = int(sys.argv[5]) if len(sys.argv) > 5 else 4_000_000
KF = float(sys.argv[6]) if len(sys.argv) > 6 else 12
LNSR = int(sys.argv[7]) if len(sys.argv) > 7 else 4
LNST = float(sys.argv[8]) if len(sys.argv) > 8 else 120
logf = open(f"grow_{tag}.log", "a")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); logf.write(s + "\n"); logf.flush()


init_col, tri = None, None
if seed == "H":
    P0 = np.array(list(heptagon().values())); P0 = unique_rows(np.concatenate([P0, conj_rows(P0)]))
elif seed == "HP":
    H = heptagon(); P0 = np.array(list(H.values())) - H["P0"]; P0 = unique_rows(np.concatenate([P0, conj_rows(P0)]))
elif seed == "G1":
    G, A, B = g1_points(); P0 = unique_rows(np.concatenate([G, conj_rows(G)]))
elif seed == "B2":
    B = sums_ball(D, 2); P0 = unique_rows(np.concatenate([B, conj_rows(B)]))
elif seed == "B1U":
    P0 = sums_ball(U, 1)
elif seed.endswith(".npy"):
    P0 = np.load(seed)
    stf = seed.replace("_pts.npy", "_state.json")
    if os.path.exists(stf):
        st = json.load(open(stf))
        tri = tuple(st.get("triangle")) if st.get("triangle") else None
        init_col = np.array(st["col"]) if st.get("col") else None
else:
    raise SystemExit("unknown seed")


def tnorm(A):
    V = (np.asarray(A, dtype=float) @ EMB) / 7.0
    return (np.abs(V) ** 2).sum(axis=1)


P = np.array(P0, dtype=np.int64)
E, _ = build_edges(P, U)
if tri is None:
    S = PointSet(P)
    idx = [int(S.index([t])[0]) for t in ([0] * 12, [7] + [0] * 11, list(D[14]))]
    tri = tuple(idx) if min(idx) >= 0 else None
    if tri is None:
        adj = [set() for _ in range(len(P))]
        for a, b in E: adj[a].add(int(b)); adj[b].add(int(a))
        for a in np.argsort(tnorm(P)):
            for b in sorted(adj[a]):
                com = adj[a] & adj[b]
                if com: tri = (int(a), int(b), int(min(com))); break
            if tri: break
tri = tuple(int(x) for x in tri)
log(f"=== grow4 seed={seed} tag={tag}: {len(P)} vertices, {len(E)} edges; triangle {tri}; maxn {MAXN}; kissat {TL}s; tabu {TABU}")


def kissat_colour(P, E, tri, tl):
    n = len(P)
    path = f"/dev/shm/grow3_{tag}.cnf"
    cls = cnf_clauses(n, E.tolist(), tri)
    with open(path, "w") as f:
        f.write(f"p cnf {4 * n} {len(cls)}\n")
        f.write("".join(" ".join(map(str, c)) + " 0\n" for c in cls))
    t = time.time()
    out = subprocess.run(["nice", "-n", "19", KISSAT, f"--time={tl}", path], capture_output=True, text=True).stdout
    os.unlink(path)
    lines = out.splitlines()
    st = next((l for l in lines if l.startswith("s ")), "s UNKNOWN")
    if "UNSATISFIABLE" in st:
        return False, None, time.time() - t
    if "SATISFIABLE" in st:
        lits = set(int(x) for l in lines if l.startswith("v ") for x in l.split()[1:] if int(x) > 0)
        col = np.array([next(c for c in range(4) if 1 + 4 * v + c in lits) for v in range(n)], dtype=np.int64)
        assert np.all(col[E[:, 0]] != col[E[:, 1]])
        return True, col, time.time() - t
    return None, None, time.time() - t


def save(status, col):
    np.save(f"grow_{tag}_pts.npy", P)
    json.dump({"seed": seed, "n": len(P), "edges": len(E), "triangle": list(tri), "status": status,
               "col": None if col is None else [int(c) for c in col]}, open(f"grow_{tag}_state.json", "w"))


col = init_col
rnd, T0 = 0, time.time()
while True:
    rnd += 1
    n = len(P)
    t1 = time.time()
    init = -np.ones(n, dtype=np.int64)
    if col is not None:
        init[:len(col)] = col
    c2 = None
    how = "lns"
    if col is not None and len(col) < n:
        adjl = [set() for _ in range(n)]
        for a, b in E:
            adjl[a].add(int(b)); adjl[b].add(int(a))
        out = lns_repair(n, adjl, list(col), list(range(len(col), n)), tri, rmax=LNSR, tl=LNST, log=lambda *a: None)
        if out is not None:
            c2 = np.array(out, dtype=np.int64)
            assert np.all(c2[E[:, 0]] != c2[E[:, 1]])
    if c2 is None:
        how = "tabu"
        c2 = tabucol(n, E, init=init, fixed=[(v, i) for i, v in enumerate(tri)], maxiter=TABU, seed=rnd)
    if c2 is None:
        how = "kissat"
        r, c2, _ = kissat_colour(P, E, tri, TL)
        if r is False:
            log(f"round {rnd}: n={n} e={len(E)}: kissat UNSAT in {time.time() - t1:.1f}s  [total {time.time() - T0:.0f}s]")
            save("unsat", None); break
        if r is None:
            log(f"round {rnd}: n={n} e={len(E)}: kissat TIME-OUT ({TL}s)  [total {time.time() - T0:.0f}s]")
            save("timeout", col); break
    col = c2
    dt = time.time() - t1
    if n >= MAXN:
        log(f"round {rnd}: n={n}: 4-colourable, reached maxn"); save("maxn", col); break
    S = PointSet(P)
    allk, allv = [], []
    for j in range(len(Ua)):
        W = P + Ua[j]
        out = np.nonzero(S.index(W) < 0)[0]
        allk.append(keys(W[out])); allv.append(np.stack([out, np.full(len(out), j)], 1))
    allk = np.concatenate(allk); allv = np.concatenate(allv)
    uk, inv, cnt = np.unique(allk, return_inverse=True, return_counts=True)
    mask = np.zeros(len(uk), dtype=np.int64)
    np.bitwise_or.at(mask, inv, 1 << col[allv[:, 0]])
    rep = np.zeros(len(uk), dtype=np.int64); rep[inv] = np.arange(len(inv))
    cand = P[allv[rep, 0]] + Ua[allv[rep, 1]]
    blocked = np.nonzero(mask == 15)[0]
    ncol = np.array([bin(m).count("1") for m in range(16)])[mask]
    if len(blocked) == 0:
        pool = np.nonzero(cnt >= 2)[0]
        if len(pool) == 0: pool = np.arange(len(uk))
        keyv = [-ncol[pool], -cnt[pool]]
    else:
        pool = blocked; keyv = [-cnt[pool]]
    inO = np.all(cand[pool] % 7 == 0, axis=1).astype(int)
    order = np.lexsort(tuple([tnorm(cand[pool]), -inO] + keyv[::-1]))
    K = max(30, min(400, int(n / KF)))
    chosen = pool[order[:K]]
    newpts = cand[chosen]
    P = np.concatenate([P, newpts])
    E, _ = build_edges(P, U)
    log(f"round {rnd}: n={n} SAT ({how}) in {dt:.1f}s; cand {len(uk)}, blocked {len(blocked)}; added {len(chosen)} "
        f"(nbrs {cnt[chosen].min()}-{cnt[chosen].max()}, inO {int(np.all(newpts % 7 == 0, axis=1).sum())}) -> n={len(P)} e={len(E)}  [total {time.time() - T0:.0f}s]")
    save("running", col)
