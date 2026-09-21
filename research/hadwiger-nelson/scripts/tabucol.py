"""TabuCol, written properly this time, and calibrated before it is believed.

The first local search here picked a RANDOM conflicting vertex and moved it to
its best colour.  That is not TabuCol and it does not work: on "Sa at five with
class 4/9 forbidden" -- which cadical proves satisfiable in ten seconds -- it
left seventy violations.  A search that cannot solve a known-satisfiable
instance says nothing about an unknown one, so its readings on the harder
instances were discarded rather than reported.

TabuCol evaluates EVERY (conflicting vertex, colour) move each step and takes
the best non-tabu one, with tenure proportional to the current conflict count.
That is the difference, and it is the whole difference.

Calibration first, on three instances whose answers are known:
  Sa at five, nothing forbidden      -- trivially satisfiable
  Sa at five, 4/9 forbidden          -- satisfiable, cadical took 10s
  Sa at four, 4/9 forbidden          -- UNSATISFIABLE, cadical took 8s
A search that solves the first two and plateaus on the third is calibrated.
Only then are the open instances worth running through it.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis

t0 = time.time()


def instance(P, classes):
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
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) in classes:
                extra.append((i, j))
    return g.n, np.array(sorted(E) + extra, dtype=np.int64)


def tabucol(n, edges, k, seed=0, iters=200_000):
    rng = np.random.default_rng(seed)
    A = np.zeros((n, n), dtype=np.int8)
    A[edges[:, 0], edges[:, 1]] = 1
    A[edges[:, 1], edges[:, 0]] = 1
    col = rng.integers(0, k, n)
    onehot = np.zeros((n, k), dtype=np.int64)
    onehot[np.arange(n), col] = 1
    gam = A @ onehot                      # gam[v, c] = neighbours of v at c
    idx = np.arange(n)
    f = int(gam[idx, col].sum()) // 2
    best = f
    tab = np.zeros((n, k), dtype=np.int64)
    for it in range(iters):
        if f == 0:
            return 0, it
        V = np.nonzero(gam[idx, col])[0]
        delta = gam[V] - gam[V, col[V]][:, None]
        allowed = (tab[V] <= it) | (f + delta < best)
        allowed[np.arange(len(V)), col[V]] = False
        if not allowed.any():
            tab[:] = 0
            continue
        cost = np.where(allowed, delta, 1 << 40) + rng.random(delta.shape)
        flat = int(np.argmin(cost))
        r, c = divmod(flat, k)
        v = int(V[r])
        cur = int(col[v])
        f += int(delta[r, c])
        tab[v, cur] = it + int(0.6 * len(V)) + int(rng.integers(0, 10)) + 1
        col[v] = c
        gam[:, cur] -= A[:, v]
        gam[:, c] += A[:, v]
        if f < best:
            best = f
    return best, iters


def run(label, P, classes, k, seeds=3, expect=None):
    n, E = instance(P, classes)
    got = []
    for sd in range(seeds):
        b, it = tabucol(n, E, k, seed=sd)
        got.append(b)
        if b == 0:
            break
    lo = min(got)
    verdict = ("found a colouring -> SATISFIABLE" if lo == 0
               else f"no colouring, best {lo} conflicts")
    flag = ""
    if expect == "SAT":
        flag = "   [calibration: OK]" if lo == 0 else "   [calibration: FAILED]"
    if expect == "UNSAT":
        flag = "   [calibration: OK]" if lo > 0 else "   [calibration: BROKEN]"
    print(f"{label}: {n} pts, {len(E)} constraints, k={k} -> {verdict}{flag}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    return lo


Sa, Y = build_Sa(F), build_Y(F)
print("-- calibration --", flush=True)
run("Sa, nothing forbidden, k=5", Sa, set(), 5, expect="SAT")
run("Sa, 4/9 forbidden, k=5", Sa, {Fr(4, 9)}, 5, expect="SAT")
run("Sa, 4/9 forbidden, k=4", Sa, {Fr(4, 9)}, 4, expect="UNSAT")
print("\n-- the open instances --", flush=True)
FOUR = {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}
run("Sa, all four classes, k=5", Sa, FOUR, 5)
run("Sa, {4/9,16/9}, k=5", Sa, {Fr(4, 9), Fr(16, 9)}, 5)
run("Y, all four classes, k=5", Y, FOUR, 5)
run("Y, 4/9, k=5", Y, {Fr(4, 9)}, 5)
run("G, 4/9, k=5", build_G(F, as_graph=False), {Fr(4, 9)}, 5)
print("DONE", flush=True)
