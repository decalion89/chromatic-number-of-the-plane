"""From a circular colouring to an orientation, and then downhill.

A (p,q)-colouring hands over a good orientation for nothing.  Order the
vertices by their position and orient each edge upward.  Around any cycle
the linear differences sum to zero, every edge has |difference| between 1
and r-1 in the scaled colouring, so the larger side of the cycle is at most
(r-1) times the smaller and |C| / min is at most r.  The orientation is
acyclic because it comes from a linear order, and no two vertices sharing a
position are adjacent, so ties are free.

That seeds the search at the value the solver already proved and lets local
search ask the next question directly: can any orientation of G score below
4.4?  If one can, chi_c(G) <= 22/5 and the pending UNSAT run is going to
come back SAT.  The step is guided rather than random -- find a cycle
witnessing the current score and flip one of its edges, since only edges on
a critical cycle can lower it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
CAR = sys.argv[1]
SEED_R = Fr(*map(int, sys.argv[2].split("/")))
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 300
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, m = len(P), len(E)
EA = np.array([a for a, _ in E])
EB = np.array([c for _, c in E])
p, q = SEED_R.numerator, SEED_R.denominator
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
sv = Solver(name="cd15", bootstrap_with=cls)
assert sv.solve(), f"{CAR} does not map to K({p}/{q})"
mod = sv.get_model()
pos = np.zeros(n, dtype=np.int64)
for v in range(n):
    for j in range(p):
        if mod[v * p + j] > 0:
            pos[v] = j
            break
print(f"{CAR}: {n} pts, {m} edges; K({p}/{q}) map found, positions used "
      f"{len(set(pos.tolist()))}/{p}  [{time.time()-t0:.0f}s]", flush=True)
orient = pos[EA] < pos[EB]


def arcs(orient, r):
    src = np.where(orient, EA, EB)
    dst = np.where(orient, EB, EA)
    u = np.concatenate([src, dst])
    v = np.concatenate([dst, src])
    w = np.concatenate([np.full(m, r - 1.0), np.full(m, -1.0)])
    return u, v, w


def neg_cycle(orient, r, want_cycle=False):
    u, v, w = arcs(orient, r)
    d = np.zeros(n)
    changed = True
    for _ in range(n + 2):
        nd = d.copy()
        np.minimum.at(nd, v, d[u] + w)
        changed = not np.array_equal(nd, d)
        d = nd
        if not changed:
            break
    if not changed:
        return None if want_cycle else False
    if not want_cycle:
        return True
    ok = d[u] + w <= d[v] + 1e-9          # admissible arcs
    nxt = {}
    for i in np.nonzero(ok)[0]:
        nxt.setdefault(int(u[i]), []).append((int(v[i]), int(i)))
    colour, stack, out = {}, [], None
    for s0 in list(nxt):
        if colour.get(s0):
            continue
        stack = [(s0, iter(nxt.get(s0, [])))]
        colour[s0], path = 1, []
        while stack:
            x, it = stack[-1]
            nx = next(it, None)
            if nx is None:
                colour[x] = 2
                stack.pop()
                if path:
                    path.pop()
                continue
            y, ai = nx
            path.append(ai)
            if colour.get(y) == 1:
                k = len(path) - 1
                cyc = [path[k]]
                while True:
                    a0 = cyc[-1]
                    if int(u[a0]) == y:
                        break
                    k -= 1
                    if k < 0:
                        break
                    cyc.append(path[k])
                return [int(a % m) for a in cyc]
            if colour.get(y) is None:
                colour[y] = 1
                stack.append((y, iter(nxt.get(y, []))))
            else:
                path.pop()
    return out


def score(orient, lo=1.0, hi=None, tol=1e-3):
    hi = float(n) if hi is None else hi
    if neg_cycle(orient, hi):
        return float("inf")
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if neg_cycle(orient, mid):
            lo = mid
        else:
            hi = mid
    return hi


cur = score(orient)
print(f"seed orientation scores {cur:.4f}  (theory says <= {float(SEED_R)})"
      f"  [{time.time()-t0:.0f}s]", flush=True)
random.seed(5)
best = cur
for step in range(STEPS):
    cyc = neg_cycle(orient, cur - 2e-3, want_cycle=True)
    if not cyc:
        print(f"   no critical cycle recovered at step {step}", flush=True)
        break
    e = random.choice(cyc)
    trial = orient.copy()
    trial[e] = not trial[e]
    s = score(trial)
    if s <= cur + 1e-9:
        orient, cur = trial, s
        if s < best - 1e-4:
            best = s
            print(f"   step {step}: {s:.4f}  (cycle had {len(cyc)} edges)"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nbest: chi_c({CAR}) <= {best:.4f}  [{time.time()-t0:.0f}s]",
      flush=True)
for r in (Fr(22, 5), Fr(13, 3), Fr(9, 2)):
    print(f"   vs K({r}) = {float(r):.4f}: "
          f"{'BELOW, so it maps' if best <= float(r) + 1e-3 else 'not reached'}",
          flush=True)
print("DONE", flush=True)
