"""Localise the obstruction, and look for a small witness.

Sixty runs of the calibrated TabuCol on Sa at five with classes 4/9 and 16/9
forbidden: none reached zero, and ONE conflict was reached eighteen times.  A
satisfiable instance that a search touches at distance one that often is
usually solved; this one never is.  Three CDCL solvers are still on it.

Two ways to get a proof faster than waiting.

BALLS.  Forbidding is monotone, so a ball about the origin that fails to
colour proves the whole thing fails and is a far smaller object to prove it
on.  Sa is a dihedral closure about the origin, so balls there respect its
symmetry.

LOCALISE.  Every near-colouring violates exactly one constraint.  Which one?
If the violated pair is the same handful across many runs, the obstruction is
local, and the subgraph induced on those vertices and their neighbourhoods is
small enough for CDCL to settle in seconds rather than hours.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

k = 5
CS = {Fr(4, 9), Fr(16, 9)}
t0 = time.time()
P0 = build_Sa(F)


def build(P):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    extra = []
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) in CS:
                extra.append((i, j))
    return g.n, sorted(E), extra


def sat(n, allpairs, budget=None):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in allpairs:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


# --- balls ----------------------------------------------------------------
O = Point(F.zero(), F.zero())
r2 = np.array([float(p.dist2(O)) for p in P0])
order = np.argsort(r2)
print("balls about the origin, classes {4/9, 16/9} forbidden:", flush=True)
for m in (100, 150, 200, 240, 280, 310, 340, 370, len(P0)):
    if m > len(P0):
        break
    Q = [P0[int(i)] for i in sorted(order[:m])]
    n, E, X = build(Q)
    t1 = time.time()
    ok = sat(n, E + X)
    print(f"  ball of {n}: {len(E)} unit + {len(X)} forbidden -> "
          f"{'colours' if ok else '*** DOES NOT COLOUR -- WITNESS ***'} in "
          f"{time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        break

# --- localise -------------------------------------------------------------
from tabucol import tabucol
n, E, X = build(P0)
edges = np.array(E + X, dtype=np.int64)
print(f"\nwhere the single conflict sits, over many runs "
      f"({len(E)} unit + {len(X)} forbidden):", flush=True)
hits = Counter()
for sd in range(40):
    A = np.zeros((n, n), dtype=np.int8)
    A[edges[:, 0], edges[:, 1]] = 1
    A[edges[:, 1], edges[:, 0]] = 1
    b, it, col = None, None, None
    rng = np.random.default_rng(5000 + sd)
    # re-run tabucol but keep the assignment at its best point
    colr = rng.integers(0, k, n)
    onehot = np.zeros((n, k), dtype=np.int64)
    onehot[np.arange(n), colr] = 1
    gam = A @ onehot
    idx = np.arange(n)
    f = int(gam[idx, colr].sum()) // 2
    best, bestcol = f, colr.copy()
    tab = np.zeros((n, k), dtype=np.int64)
    for it in range(300_000):
        if f == 0:
            break
        V = np.nonzero(gam[idx, colr])[0]
        delta = gam[V] - gam[V, colr[V]][:, None]
        allowed = (tab[V] <= it) | (f + delta < best)
        allowed[np.arange(len(V)), colr[V]] = False
        if not allowed.any():
            tab[:] = 0
            continue
        cost = np.where(allowed, delta, 1 << 40) + rng.random(delta.shape)
        r, c = divmod(int(np.argmin(cost)), k)
        v = int(V[r])
        cur = int(colr[v])
        f += int(delta[r, c])
        tab[v, cur] = it + int(0.6 * len(V)) + int(rng.integers(0, 10)) + 1
        colr[v] = c
        gam[:, cur] -= A[:, v]
        gam[:, c] += A[:, v]
        if f < best:
            best, bestcol = f, colr.copy()
    if best == 1:
        bad = [(a, b) for a, b in E + X if bestcol[a] == bestcol[b]]
        for e in bad:
            hits[e] += 1
    if sd % 10 == 9:
        print(f"  {sd+1} runs, {len(hits)} distinct violated pairs, "
              f"top {hits.most_common(6)}  [{time.time()-t0:.0f}s]",
              flush=True)
print(f"\n{len(hits)} distinct pairs ever left violated; "
      f"top {hits.most_common(10)}", flush=True)
print("DONE", flush=True)
