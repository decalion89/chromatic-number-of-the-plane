"""min3multi.py -- min3inc.py repeated with random deletion orders (same incremental solver), keeping the
smallest vertex-critical graph found.

usage: python3 min3multi.py STATE.json PTS.npy OUT_PREFIX RUNS [TIME_LIMIT_S]

min3inc.py: shrink a non-3-colourable unit-distance graph to a vertex-critical one with ONE incremental
SAT solver (CaDiCaL 1.5.3 through pysat): a selector literal per vertex switches its at-least-one-colour
clause on, each test is a solve under assumptions, and an UNSAT answer returns a core of selectors, i.e.
a subset of the vertices that is still not 3-colourable. Faster than min3.py on large graphs.

usage: python3 min3inc.py STATE.json PTS.npy OUT_PREFIX
1. S = the 3-core; S = the core of the first refutation, then its 3-core.
2. for each vertex v, lowest degree first: if the 3-core of S - v is not 3-colourable, S = the core of that
   refutation (3-core); otherwise v is necessary (the solver found a 3-colouring of S - v without v's
   unused neighbours). Necessity is monotone, so at the end every vertex of S is necessary.
Writes OUT_PREFIX_pts.npy and OUT_PREFIX.json like min3.py; certify_q.py checks the result from scratch."""
import sys, json, time
import numpy as np
from pysat.solvers import Solver

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


def x(v, c):
    return 3 * v + c + 1


def sel(v):
    return 3 * n + v + 1


S = core3(range(n))
T0 = time.time()
print(f"start {n}, 3-core {len(S)}", flush=True)
solver = Solver(name="cadical153")
for v in S:
    solver.add_clause([-sel(v), x(v, 0), x(v, 1), x(v, 2)])
    for w in adj[v]:
        if w in S and v < w:
            for c in range(3):
                solver.add_clause([-x(v, c), -x(w, c)])


def refute(T):
    """None if T is 3-colourable, else the vertices of a core of the refutation"""
    if solver.solve(assumptions=[sel(v) for v in sorted(T)]):
        return None
    return {l - 3 * n - 1 for l in solver.get_core() if l > 3 * n}


C = refute(S)
assert C is not None, "the 3-core is 3-colourable"
S = core3(C)
print(f"first core: {len(S)}  [{time.time() - T0:.0f}s]", flush=True)
import random
RUNS = int(sys.argv[4]) if len(sys.argv) > 4 else 50
TLIM = float(sys.argv[5]) if len(sys.argv) > 5 else 3600
S0 = set(S)
best = None
for run in range(RUNS):
    if time.time() - T0 > TLIM:
        break
    rng = random.Random(run)
    S = set(S0)
    if run:                                     # a fresh core from a random subset order is not available; shuffle
        order = sorted(S); rng.shuffle(order)
        order.sort(key=lambda v: len(adj[v] & S) + rng.random() * 3)
    else:
        order = sorted(S, key=lambda v: (len(adj[v] & S), v))
    necessary = set()
    tests = 0
    for v in order:
        if v not in S or v in necessary:
            continue
        tests += 1
        T = core3(S - {v})
        C = refute(T)
        if C is None:
            necessary.add(v)
        else:
            S = core3(C)
    assert necessary >= S
    if best is None or len(S) < len(best):
        best = set(S)
    print(f"run {run}: {len(S)} vertices (best {len(best)})  [{time.time() - T0:.0f}s]", flush=True)
S = best
print(f"vertex-critical: {len(S)} vertices", flush=True)
L = sorted(S)
pos = {v: k for k, v in enumerate(L)}
E = [[pos[a], pos[b]] for a in L for b in adj[a] if b in pos and a < b]
np.save(out + "_pts.npy", np.array([pts[v] for v in L], dtype=np.int64))
json.dump({"D": D, "U": [list(u) for u in U], "points": [list(pts[v]) for v in L], "edges": E}, open(out + ".json", "w"))
print(f"saved {out}: {len(L)} vertices, {len(E)} edges", flush=True)
