"""chi_c from above, by orientations rather than by SAT.

For a fixed orientation the maximum of |C| / |C+| over cycles is a
min-ratio-cycle problem: give an edge traversed with the orientation the
weight r-1 and one traversed against it the weight -1, so a cycle totals
r|C+| - |C| and is negative exactly when its ratio exceeds r.  Bellman-Ford
finds such a cycle in polynomial time, and binary search on r gives the
orientation's exact score.  Any orientation's score is an upper bound on
chi_c, and it is a certificate: no search has to be exhaustive for the bound
to hold.

That is the opposite trade from the circular-clique encoding, which proves
lower bounds by UNSAT and chokes on size.  Here the graph can be large and
the answer comes in seconds, but it can only ever say "no higher than".
Used to eliminate candidates, and to settle quickly whether a pending UNSAT
run is in fact going to come back SAT.
"""
import sys, time, pickle, random
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1]
TRIES = int(sys.argv[2]) if len(sys.argv) > 2 else 40
t0 = time.time()
if CAR in ("G", "Sa", "Y"):
    P = {"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
         "Y": build_Y}[CAR](K)
else:
    P = pickle.load(open(SC + CAR, "rb"))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, m = len(P), len(E)
EA = np.array([a for a, _ in E], dtype=np.int64)
EB = np.array([c for _, c in E], dtype=np.int64)
print(f"{CAR}: {n} points, {m} edges  [{time.time()-t0:.0f}s]", flush=True)


def has_ratio_above(orient, r):
    """Is there a cycle with |C| / |C+| > r?  Bellman-Ford, vectorised."""
    src = np.where(orient, EA, EB)
    dst = np.where(orient, EB, EA)
    # with the orientation: weight r-1 ; against it: weight -1
    u = np.concatenate([src, dst])
    v = np.concatenate([dst, src])
    w = np.concatenate([np.full(m, r - 1.0), np.full(m, -1.0)])
    d = np.zeros(n)
    for _ in range(n):
        cand = d[u] + w
        nd = d.copy()
        np.minimum.at(nd, v, cand)
        if np.allclose(nd, d, atol=1e-12):
            return False
        d = nd
    return True


def score(orient, lo=1.0, hi=None, tol=1e-3):
    hi = float(n) if hi is None else hi
    if not has_ratio_above(orient, hi):
        while hi - lo > tol:
            mid = (lo + hi) / 2
            if has_ratio_above(orient, mid):
                lo = mid
            else:
                hi = mid
        return hi
    return float("inf")


def from_order(perm):
    """Acyclic orientation induced by a vertex ordering."""
    pos = np.empty(n, dtype=np.int64)
    pos[perm] = np.arange(n)
    return pos[EA] < pos[EB]


# A random vertex order is a terrible orientation: on hundreds of points it
# leaves long cycles with a single minority edge, so its ratio is the cycle
# length.  Proper colourings give good ones for free -- orienting from lower
# class to higher makes every cycle's minority side at least |C|/k, which is
# Minty's direction of the theorem -- so the search starts there and shuffles
# the class order, which changes the orientation without changing the bound.
from pysat.solvers import Solver
KC = int(sys.argv[3]) if len(sys.argv) > 3 else 5
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, c in E:
    for col in range(KC):
        cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
sv = Solver(name="cd15", bootstrap_with=cls)
assert sv.solve(), f"{CAR} is not {KC}-colourable"
mod = sv.get_model()
colour = np.zeros(n, dtype=np.int64)
for v in range(n):
    for c in range(KC):
        if mod[v * KC + c] > 0:
            colour[v] = c
            break
print(f"{KC}-colouring found, classes "
      f"{[int((colour == c).sum()) for c in range(KC)]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
random.seed(3)
best = float("inf")
for t in range(TRIES):
    if t == 0:
        rank = np.arange(KC)
    else:
        rank = np.array(random.sample(range(KC), KC))
    key = rank[colour] * (n + 1) + np.array(
        random.sample(range(n), n) if t else np.arange(n))
    orient = key[EA] < key[EB]
    sc = score(orient)
    if sc < best:
        best = sc
        print(f"   try {t}: {sc:.4f}   (chi_c <= {sc:.4f})"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nbest orientation found: chi_c({CAR}) <= {best:.4f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for r in (Fr(22, 5), Fr(9, 2), Fr(5)):
    print(f"   vs K({r}) = {float(r):.4f}: "
          f"{'bound is below, so it MAPS' if best <= float(r) + 1e-3 else 'no conclusion from above'}",
          flush=True)
print("DONE", flush=True)
