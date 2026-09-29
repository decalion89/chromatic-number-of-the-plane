"""grow.py -- colouring-guided growth of a unit-distance graph with edges in U (126 directions).

usage: python3 grow.py SEED TAG [maxn] [tlimit_per_solve]
  SEED: H (heptagon + mirror), HP (heptagon moved to P0=0, + mirror), G1 (G1 u conj G1), B2 (B_2(D) u conj),
        B1U (B_1(U)), or a .npy file of points
Every round: take the current 4-colouring c; candidates y = x + u (x in graph, u in U) not in the graph;
"blocked" candidates see all 4 colours among their graph-neighbours; add the best of them
(most neighbours, then in Z[zeta21], then smallest trace norm), re-solve incrementally (CaDiCaL 1.9.5,
phases from c). Stops at UNSAT, at maxn vertices, or at a solver time-out. Saves the state every round.
"""
import sys, time, json
import numpy as np
from flat import *
from seeds import heptagon, g1_points, EMB
from satutil import solve_limited, write_dimacs
from tabu import tabucol
from pysat.solvers import Solver

D, CD, U = directions()
Ua = np.array(U, dtype=np.int64)
seed, tag = sys.argv[1], sys.argv[2]
MAXN = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
TL = float(sys.argv[4]) if len(sys.argv) > 4 else 900
TABU = int(sys.argv[5]) if len(sys.argv) > 5 else 4_000_000
logf = open(f"grow_{tag}.log", "a")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    logf.write(s + "\n"); logf.flush()


if seed == "H":
    P0 = np.array(list(heptagon().values()))
    P0 = unique_rows(np.concatenate([P0, conj_rows(P0)]))
elif seed == "HP":
    H = heptagon(); P0 = np.array(list(H.values())) - H["P0"]
    P0 = unique_rows(np.concatenate([P0, conj_rows(P0)]))
elif seed == "G1":
    G, A, B = g1_points()
    P0 = unique_rows(np.concatenate([G, conj_rows(G)]))
elif seed == "B2":
    B = sums_ball(D, 2)
    P0 = unique_rows(np.concatenate([B, conj_rows(B)]))
elif seed == "B1U":
    P0 = sums_ball(U, 1)
elif seed.endswith(".npy"):
    P0 = np.load(seed)
else:
    raise SystemExit("unknown seed")


def tnorm(A):
    """trace norm sum_sigma |sigma(x)|^2 of x = A/7"""
    V = (np.asarray(A, dtype=float) @ EMB) / 7.0
    return (np.abs(V) ** 2).sum(axis=1)


class G:
    pass


g = G()
g.P = []              # list of np arrays
g.idx = {}            # tuple -> index
g.adj = []
g.s = Solver(name="cadical195")
g.nedges = 0
g.E = []
g.col = None


def var(v, c):
    return 1 + 4 * v + c


def add_points(pts):
    new = []
    for p in pts:
        t = tuple(int(x) for x in p)
        if t in g.idx:
            continue
        v = len(g.P)
        g.idx[t] = v
        g.P.append(np.array(t, dtype=np.int64))
        g.adj.append([])
        g.s.add_clause([var(v, c) for c in range(4)])
        new.append(v)
    for v in new:
        p = g.P[v]
        W = p + Ua
        for w in W:
            j = g.idx.get(tuple(int(x) for x in w))
            if j is not None and j != v and j not in g.adj[v]:
                g.adj[v].append(j); g.adj[j].append(v)
                g.nedges += 1
                g.E.append((v, j))
                for c in range(4):
                    g.s.add_clause([-var(v, c), -var(j, c)])
    return new


add_points(P0)
# fix a triangle: prefer 0,1,zeta6 if present, else any triangle containing a vertex nearest the origin
tri = None
t0p = (tuple([0] * 12), tuple([7] + [0] * 11), D[14])
if all(t in g.idx for t in t0p):
    tri = tuple(g.idx[t] for t in t0p)
else:
    order = np.argsort(tnorm(np.array(g.P)))
    for a in order:
        for b in g.adj[a]:
            com = set(g.adj[a]) & set(g.adj[b])
            if com:
                tri = (int(a), int(b), int(min(com))); break
        if tri: break
assert tri is not None
for i, v in enumerate(tri):
    g.s.add_clause([var(v, i)])
log(f"=== grow seed={seed} tag={tag}: {len(g.P)} vertices, {g.nedges} edges; fixed triangle {tri}; maxn {MAXN}, tlimit {TL}s")


def save(status):
    np.save(f"grow_{tag}_pts.npy", np.array(g.P))
    json.dump({"seed": seed, "n": len(g.P), "edges": g.nedges, "triangle": tri, "status": status},
              open(f"grow_{tag}_state.json", "w"))


rnd = 0
T0 = time.time()
while True:
    rnd += 1
    n = len(g.P)
    t1 = time.time()
    init = -np.ones(n, dtype=np.int64)
    if g.col is not None:
        init[:len(g.col)] = g.col
    colt = tabucol(n, np.array(g.E, dtype=np.int64), init=init, fixed=[(v, i) for i, v in enumerate(tri)],
                   maxiter=TABU, seed=rnd)
    how = "tabu"
    if colt is not None:
        r = True
    else:
        how = "cadical"
        r = solve_limited(g.s, tlimit=TL)
    dt = time.time() - t1
    if r is None:
        log(f"round {rnd}: n={n} e={g.nedges}: solver TIME-OUT after {dt:.0f}s")
        save("timeout"); break
    if r is False:
        log(f"round {rnd}: n={n} e={g.nedges}: UNSAT (pysat CaDiCaL) in {dt:.1f}s  [total {time.time() - T0:.0f}s]")
        save("unsat"); break
    if colt is not None:
        col = colt
    else:
        model = g.s.get_model()
        pos = set(l for l in model if l > 0)
        col = np.array([next(c for c in range(4) if var(v, c) in pos) for v in range(n)], dtype=np.int64)
    g.col = col
    # phases: keep this colouring as preferred
    g.s.set_phases([var(v, int(col[v])) for v in range(n)])
    if n >= MAXN:
        log(f"round {rnd}: n={n}: SAT, reached maxn; stop"); save("maxn"); break
    # candidates
    P = np.array(g.P)
    S = PointSet(P)
    allk, allv = [], []
    for j in range(len(Ua)):
        W = P + Ua[j]
        idx = S.index(W)
        out = np.nonzero(idx < 0)[0]
        allk.append(keys(W[out])); allv.append(np.stack([out, np.full(len(out), j)], 1))
    allk = np.concatenate(allk); allv = np.concatenate(allv)
    uk, inv, cnt = np.unique(allk, return_inverse=True, return_counts=True)
    mask = np.zeros(len(uk), dtype=np.int64)
    np.bitwise_or.at(mask, inv, 1 << col[allv[:, 0]])
    rep = np.zeros(len(uk), dtype=np.int64); rep[inv] = np.arange(len(inv))
    reps = allv[rep]
    cand = P[reps[:, 0]] + Ua[reps[:, 1]]
    blocked = np.nonzero(mask == 15)[0]
    ncol = np.array([bin(m).count("1") for m in range(16)])[mask]
    if len(blocked) == 0:
        # nothing blocked: add the candidates with the most neighbours (3 colours seen first)
        pool = np.nonzero(cnt >= 2)[0]
        if len(pool) == 0:
            pool = np.arange(len(uk))
        keyv = (-ncol[pool], -cnt[pool])
    else:
        pool = blocked
        keyv = (-cnt[pool],)
    inO = np.all(cand[pool] % 7 == 0, axis=1).astype(int)
    tn = tnorm(cand[pool])
    # np.lexsort sorts by the last key first
    order = np.lexsort(tuple([tn, -inO] + list(keyv)[::-1]))
    K = max(30, min(400, n // 12))
    chosen = pool[order[:K]]
    new = add_points(cand[chosen])
    log(f"round {rnd}: n={n} e={g.nedges - 0} SAT ({how}) in {dt:.1f}s; candidates {len(uk)}, blocked {len(blocked)}; "
        f"added {len(new)} (nbrs {cnt[chosen].min() if len(chosen) else 0}-{cnt[chosen].max() if len(chosen) else 0}, "
        f"inO {int(np.all(cand[chosen] % 7 == 0, axis=1).sum())}) -> n={len(g.P)}  [total {time.time() - T0:.0f}s]")
    if rnd % 5 == 0:
        save("running")
