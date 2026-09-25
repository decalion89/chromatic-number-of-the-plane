"""alpha(H), and what it says about which rung H could possibly refuse.

Mapping to K(p/q) needs p independent sets covering every vertex exactly q
times, so it needs alpha(H) >= (q/p)*n.  Turn that round: a graph whose
independence ratio is below q/p refuses K(p/q) outright, no solver needed.
Measuring alpha therefore says which rungs are settled for free and which
need the structure argument -- and n/alpha is a lower bound on chi_f, which
places the graph against the fractional scale Polymath16 pushed.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1]
t0 = time.time()
P = (build_G(K, as_graph=False) if CAR == "G"
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
adj = [set() for _ in range(n)]
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
print(f"{CAR}: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]",
      flush=True)
# greedy lower bound, least-degree-first, repeated with removal
best = []
order = sorted(range(n), key=lambda v: len(adj[v]))
banned, cur = set(), []
for v in order:
    if v not in banned:
        cur.append(v)
        banned |= adj[v] | {v}
best = cur
print(f"greedy independent set: {len(best)}  ratio "
      f"{len(best)/n:.4f}  [{time.time()-t0:.0f}s]", flush=True)


def has(k):
    pool = IDPool(start_from=n + 1)
    cls = [[-(a + 1), -(c + 1)] for a, c in E]
    cls += list(CardEnc.atleast(lits=list(range(1, n + 1)), bound=k,
                                vpool=pool,
                                encoding=EncType.seqcounter).clauses)
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    return ok


lo, hi = len(best), n
while lo < hi:
    mid = (lo + hi + 1) // 2
    ok = has(mid)
    print(f"   alpha >= {mid}? {'yes' if ok else 'no'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if ok:
        lo = mid
    else:
        hi = mid - 1
print(f"\nalpha({CAR}) = {lo}, ratio {lo/n:.4f}, "
      f"chi_f >= n/alpha = {n/lo:.4f}  [{time.time()-t0:.0f}s]", flush=True)
for r in (Fr(9, 2), Fr(13, 3), Fr(5, 1)):
    need = r.denominator / r.numerator
    print(f"   K({r}) needs ratio >= {need:.4f}: "
          f"{'possible' if lo / n >= need else 'REFUSED by counting alone'}",
          flush=True)
print("DONE", flush=True)
