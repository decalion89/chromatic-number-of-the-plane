"""Shrink G to a 5-critical core, because size is what the cap is short of.

The cap at four colours is abundant because the smallest 4-critical
unit-distance graph is the Moser spindle, seven vertices, and Sa carries 228
of them over 397 points -- 0.57 per point.  Several pass through every ring
point at once, and that overlap is what holds a hexagon to two colours.

At five the smallest known 5-chromatic unit-distance graph has about 500
vertices.  Gp carries twelve copies of G over 13873 points: 0.0009 per point,
six hundred times thinner.  Nothing overlaps, nothing is rigid, and no ring is
capped.  That single ratio explains every negative in this work.

So the lever is size.  A smaller 5-chromatic graph makes a denser carrier
possible, and G at 1581 vertices is nowhere near critical -- de Grey built it
for provability, not for economy.

The cheap way down is the solver's own unsatisfiable core: the 4-colouring
formula for G has no model, and the core names a subset of clauses that
already has none, which is a 5-chromatic subgraph.  Then remove vertices one
at a time, keeping only those whose removal would make the rest 4-colourable,
until every vertex is essential.  What is left is 5-critical.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
k = 4                     # we keep the graph NOT k-colourable, i.e. chi >= 5
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
print(f"G: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)


def not_colourable(keep):
    """True when the induced subgraph on `keep` needs more than k colours."""
    loc = {v: i for i, v in enumerate(sorted(keep))}
    m = len(loc)
    if m <= k:
        return False
    cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
    for a, c in E:
        if a in keep and c in keep:
            for col in range(k):
                cls.append([-(1 + loc[a] * k + col), -(1 + loc[c] * k + col)])
    for col in range(1, k):
        cls.append([-(1 + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    out = s.solve()
    s.delete()
    return not out


keep = set(range(n))
t1 = time.time()
assert not_colourable(keep), "G is 4-colourable?!"
print(f"confirmed chi(G) > 4 in {time.time()-t1:.1f}s  [{time.time()-t0:.0f}s]",
      flush=True)

# Cheap first pass: anything with degree below k cannot matter.
changed = True
while changed:
    changed = False
    d = {v: sum(1 for u in adj[v] if u in keep) for v in keep}
    low = {v for v in keep if d[v] < k}
    if low:
        keep -= low
        changed = True
print(f"after stripping degree < {k}: {len(keep)} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)
assert not_colourable(keep), "stripping broke it"

order = sorted(keep, key=lambda v: sum(1 for u in adj[v] if u in keep))
removed = 0
for v in order:
    if v not in keep:
        continue
    cand = keep - {v}
    if not_colourable(cand):
        keep = cand
        removed += 1
        if removed % 50 == 0:
            print(f"   removed {removed}, {len(keep)} left"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
kept = sorted(keep)
ke = [(a, c) for a, c in E if a in keep and c in keep]
print(f"\n5-critical core of G: {len(kept)} points, {len(ke)} edges "
      f"({2*len(ke)/max(len(kept),1):.2f}/v), down from {n}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/Gcore.pkl",
          "wb") as f:
    pickle.dump([P[i] for i in kept], f)
print("DONE", flush=True)
