"""Local search, to tell a hard SAT from an UNSAT without waiting for either.

CDCL explodes on both sides of a phase transition, so "three hours and no
answer" says the instance is hard and nothing about which way it goes.  Local
search is asymmetric: on a satisfiable instance a tabu search over colourings
usually lands on a solution in seconds, and on an unsatisfiable one it
plateaus above zero forever.  That is evidence rather than proof in the second
case, but it is the right kind of evidence and it arrives immediately.

The formulation is the natural one: a colouring is an assignment of five
colours to the vertices, its cost is the number of violated constraints --
graph edges plus the forbidden pairs, which enter on exactly the same footing
since forbidding a pair IS adding the edge -- and a move recolours one vertex.
Tabu on (vertex, colour) for a few iterations stops it cycling.

Run on the instances CDCL has not answered: Sa at five with all four carrying
classes, and with the pair {4/9, 16/9}, and G at five with 4/9 alone.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis

k = 5
t0 = time.time()


def instance(P, classes):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    pr = []
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) in classes:
                pr.append((i, j))
    return g.n, np.array(sorted(E) + pr, dtype=np.int32)


def tabu(n, edges, seed=0, iters=400_000, tenure=10):
    """Tabu over colourings: cost is the number of violated constraints."""
    rng = np.random.default_rng(seed)
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    adj = [np.array(x, dtype=np.int64) for x in adj]
    col = rng.integers(0, k, n)
    cnt = np.zeros((n, k), dtype=np.int64)
    for v in range(n):
        if adj[v].size:
            np.add.at(cnt[v], col[adj[v]], 1)
    idx = np.arange(n)
    cost = int(cnt[idx, col].sum()) // 2
    best = cost
    tab = np.zeros((n, k), dtype=np.int64)
    for it in range(iters):
        if cost == 0:
            return 0, it
        bad = np.nonzero(cnt[idx, col])[0]
        v = int(bad[rng.integers(0, bad.size)])
        cur = int(col[v])
        d = cnt[v] - cnt[v, cur]
        free = (tab[v] <= it) | (cost + d < best)
        free[cur] = False
        if not free.any():
            continue
        dd = np.where(free, d, 1 << 30)
        c = int(np.argmin(dd + rng.random(k)))
        cost += int(d[c])
        tab[v, cur] = it + tenure + int(rng.integers(0, 6))
        col[v] = c
        if adj[v].size:
            np.add.at(cnt[:, cur], adj[v], -1)
            np.add.at(cnt[:, c], adj[v], 1)
        if cost < best:
            best = cost
    return best, iters


JOBS = [
    ("Sa, all four classes", build_Sa, {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}),
    ("Sa, {4/9, 16/9}", build_Sa, {Fr(4, 9), Fr(16, 9)}),
    ("Sa, 4/9 alone", build_Sa, {Fr(4, 9)}),
    ("Y, all four classes", build_Y, {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}),
    ("G, 4/9 alone", build_G, {Fr(4, 9)}),
]
for name, build, cs in JOBS:
    P = build(F) if build is not build_G else build(F, as_graph=False)
    n, E = instance(P, cs)
    got = []
    for sd in range(3):
        b, it = tabu(n, E, seed=sd)
        got.append(b)
        if b == 0:
            break
    ok = min(got) == 0
    verdict = ("SATISFIABLE (local search found a colouring)" if ok
               else f"no colouring found, best {min(got)} violations left")
    print(f"{name}: {n} pts, {len(E)} constraints -> {verdict}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
