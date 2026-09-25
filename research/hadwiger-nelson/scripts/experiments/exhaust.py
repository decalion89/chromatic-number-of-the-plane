"""Every closable pair of G, tested exactly.  No ranking, no sampling.

The graded metric ranked pairs by how often they agreed in 32 sampled
colourings, and the null control showed that ranking is noise: two copies of G
with no cross edge at all reach 23 and 24 of 32, against 21 and 22 for G
itself.  The maximum over a million 32-sample estimates is an extreme value,
not a signal, so ranking by it ranks nothing.

Drop the ranking.  What the spindle needs is a pair that is forced to one
colour in EVERY proper 5-colouring, and that is decidable exactly: the pair
(i, j) is forced-same iff G with the edge (i, j) added is not 5-colourable,
and by colour symmetry that is the single call

    solve(assumptions=[x_{i,0}, -x_{j,0}])

returning UNSAT.  One call, no samples, no error bars.

A million calls is too many, but the spindle also needs the pair's distance to
be closable -- there is no rotation to build otherwise -- and that filter is
severe.  So: enumerate every pair of G at a RATIONAL squared distance (all
sixteen non-constant coordinates of dx^2 + dy^2 vanish), keep the closable
ones, and test each exactly.  The answer is complete for G, not sampled.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
P = build_G(F, as_graph=False)
g = build_graph(P)
n = g.n
E = set((min(a, b), max(a, b)) for a, b in g.edges())
print(f"G: {n} pts, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

basis = IntBasis.covering(P)
rows = basis.rows(P)
dim = basis.dim
D2 = basis.D * basis.D
print(f"denominator D={basis.D}, dim={dim}  [{time.time()-t0:.0f}s]", flush=True)

# --- every pair at a rational squared distance -----------------------------
pairs = []                       # (i, j, Fraction d2)
for i in range(n - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) in E:
            continue
        pairs.append((i, j, Fr(int(sq[off, 0]), D2)))
    if i % 300 == 0:
        print(f"  pivot {i}: {len(pairs)} rational so far"
              f"  [{time.time()-t0:.0f}s]", flush=True)

print(f"\n{len(pairs)} non-edge pairs at a RATIONAL squared distance"
      f"  [{time.time()-t0:.0f}s]", flush=True)
byd = Counter(p[2] for p in pairs)
print(f"{len(byd)} distinct rational squared distances", flush=True)

# --- keep the closable ones ------------------------------------------------
ok = {d: closable_distance(d) for d in byd}
clo = [p for p in pairs if ok[p[2]]]
print(f"{sum(ok.values())} of them closable, carrying {len(clo)} pairs",
      flush=True)
for d, c in sorted(byd.items(), key=lambda t: -t[1])[:20]:
    print(f"    D={d}  {c} pairs  closable={ok[d]}", flush=True)

# --- the exact test --------------------------------------------------------
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in sorted(E):
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="cd15", bootstrap_with=cls)
assert sv.solve(), "G is 5-colourable"
print(f"\ntesting {len(clo)} closable pairs exactly"
      f"  [{time.time()-t0:.0f}s]", flush=True)
forced = []
for t, (i, j, d) in enumerate(clo):
    if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)]):
        forced.append((i, j, d))
        print(f"  *** FORCED SAME: {i},{j} at D={d} ***", flush=True)
    if t and t % 2000 == 0:
        print(f"  {t}/{len(clo)}, {len(forced)} forced"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nDONE: {len(forced)} forced-same pairs among {len(clo)} closable"
      f"  [{time.time()-t0:.0f}s]", flush=True)
