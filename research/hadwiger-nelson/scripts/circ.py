"""The circular chromatic number: the same question with a scale on it.

Every measurement in this project has asked a yes-or-no question -- is this
graph 5-colourable -- and the answer has been yes about a hundred times, with
no way to tell which of two failures came closer.  That is what "there is no
gradient" has meant here all along.

There is a gradient, and it is the circular chromatic number.  A graph maps
to the circular clique K(p/q) when its vertices can be placed on p points of
a circle so that adjacent ones land at least q apart either way round, and
chi_c is the least p/q for which that works.  It satisfies chi_c <= chi and
ceil(chi_c) = chi, so it is a real number sitting inside the integer, and it
distinguishes graphs that chi cannot.

The point for this search: chi_c > 5 implies chi >= 6.  So a graph with
chi_c = 4.97 is nearer the goal than one with 4.83, and for the first time
that difference is visible.  Every 5-chromatic graph here has chi_c somewhere
in (4, 5], and nothing in this work has ever looked at where.

Encoded directly: one position variable per vertex per point of the circle,
and for each edge the forbidden window of positions around each choice.  That
is (2q-1)*p clauses per edge rather than p^2, which keeps it affordable.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pairsat import pairs_at
from pysat.solvers import Solver

SEED = sys.argv[1] if len(sys.argv) > 1 else "G"
EXTRA = sys.argv[2] if len(sys.argv) > 2 else ""
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
if SEED.endswith(".pkl"):
    P = pickle.load(open(SC + SEED, "rb"))
else:
    P = {"G": lambda: build_G(K, as_graph=False), "Sa": lambda: build_Sa(K),
         "Y": lambda: build_Y(K)}[SEED]()
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
if EXTRA:
    E |= set(pairs_at(b, r, int(Fr(EXTRA) * b.D * b.D)))
E = sorted(E)
n = len(P)
print(f"{SEED}{' + d^2=' + EXTRA if EXTRA else ''}: {n} points, {len(E)} "
      f"edges  [{time.time()-t0:.0f}s]", flush=True)


def maps_to(p, q):
    """Is there a homomorphism to the circular clique K(p/q)?"""
    cls = [[1 + v * p + i for i in range(p)] for v in range(n)]
    for v in range(n):                      # at most one position each
        for i in range(p):
            for j in range(i + 1, p):
                cls.append([-(1 + v * p + i), -(1 + v * p + j)])
    for a, c in E:
        for i in range(p):
            for d in range(-(q - 1), q):
                j = (i + d) % p
                cls.append([-(1 + a * p + i), -(1 + c * p + j)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    return ok


cands = sorted({(p, q) for q in range(1, 10) for p in range(4 * q + 1,
                                                            5 * q + 1)
                if Fr(p, q).denominator == q and 4 < Fr(p, q) <= 5},
               key=lambda pq: Fr(*pq))
cands = [pq for pq in cands if pq[0] * n <= 400000]
print(f"binary search over {len(cands)} circular cliques in (4, 5]:"
      f" {float(Fr(*cands[0])):.4f} .. {float(Fr(*cands[-1])):.4f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# The set of ratios a graph maps to is an UP-set: mapping to K(p/q) gives
# every larger ratio for free.  So the boundary is found by bisection, and
# chi_c is the least ratio that still admits a homomorphism.
lo, hi, seen = 0, len(cands) - 1, {}


def test(i):
    if i not in seen:
        p, q = cands[i]
        seen[i] = maps_to(p, q)
        print(f"   K({p}/{q}) = {float(Fr(p,q)):.4f}: "
              f"{'maps' if seen[i] else 'does NOT map'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    return seen[i]


if test(hi):
    while lo < hi:
        mid = (lo + hi) // 2
        if test(mid):
            hi = mid
        else:
            lo = mid + 1
    p, q = cands[lo]
    print(f"\nchi_c = {Fr(p, q)} = {float(Fr(p, q)):.4f}   (chi = "
          f"{-(-p // q)})", flush=True)
    print(f"distance from 5: {5 - float(Fr(p, q)):.4f}", flush=True)
else:
    print("\nno circular clique at or below 5 admits it -- chi_c > 5, "
          "so chi >= 6", flush=True)
print("DONE", flush=True)
