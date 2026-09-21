"""Cone over a circle: the construction that needs no spindle at all.

Every search here has hunted pairs forced to share a colour, because that is
what the spindle consumes.  The apex wants the opposite.  If P is a set of
points on one unit circle and no proper k-colouring of the graph squeezes P
into k-1 colours, then the circle's centre sees all k colours among its
neighbours and cannot be coloured: the cone over P needs k+1, and the cone is
still a unit-distance graph because P is at distance one from its centre.

Two measurements, both exact and both one SAT call each.

  - forced_different, the dual of forced_same: must these two vertices differ
    in every colouring?  Five such on one circle finish the job outright.
  - the squeeze: can P be confined to k-1 colours?  By colour symmetry one
    call settles it for every choice of which colour to drop.

Run over G at five, whose every vertex is a candidate centre, and over every
point of the field that G puts on a unit circle.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import forced_same, forced_different
from pysat.solvers import Solver

k = 5
t0 = time.time()
P = build_G(K, as_graph=False)
n = len(P)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
E = sorted(set((min(a, b), max(a, b)) for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
adj = defaultdict(set)
for a, b in E:
    adj[a].add(b)
    adj[b].add(a)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in E:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="cd15", bootstrap_with=cls)
assert sv.solve()
print(f"G: {n} pts, {len(E)} edges, 5-colourable  [{time.time()-t0:.0f}s]",
      flush=True)

# 1. how many ADJACENT pairs are there -- trivially forced different -- and is
#    any NON-adjacent pair forced different at all?
nonadj_tested = 0
fd = []
# test non-adjacent pairs that could share a unit circle: any pair at squared
# distance < 4 has a common point at distance one from both
cand = []
for i in range(n):
    d = rows - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = int(off)
        if j <= i:
            continue
        v = Fr(int(sq[j, 0]), D2)
        if 0 < v < 4 and (i, j) not in Eset:
            cand.append((i, j, v))
print(f"{len(cand)} non-adjacent pairs at rational squared distance < 4"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for i, j, v in cand:
    nonadj_tested += 1
    if forced_different(sv, i, j, k):
        fd.append((i, j, v))
        if len(fd) <= 5:
            print(f"*** forced DIFFERENT and non-adjacent: {i},{j} at d^2={v}"
                  f" ***", flush=True)
    if nonadj_tested % 20000 == 0:
        print(f"   {nonadj_tested}/{len(cand)}, {len(fd)} so far"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nforced-different non-adjacent pairs: {len(fd)} of {len(cand)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
if fd:
    print("   by distance:", Counter(str(v) for _, _, v in fd).most_common(10),
          flush=True)

# 2. the squeeze, over every existing vertex as the circle's centre
best = []
for c in range(n):
    Pc = sorted(adj[c])
    if len(Pc) < k:
        continue
    lits = [-(1 + p * k + 0) for p in Pc]          # nobody on the circle is 0
    ok = sv.solve(assumptions=lits)
    if not ok:
        print(f"*** CENTRE {c}: its {len(Pc)} neighbours cannot avoid a "
              f"colour -- chi >= 6 ***", flush=True)
        best.append(c)
print(f"{len(best)} centres whose neighbourhood forces all five colours"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
