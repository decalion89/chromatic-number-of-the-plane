"""The width-two hub spindle over every even number of steps, not just two.

Exhaustive transversal search corrected the hand argument: it is not only two
steps that close.  EVERY even m closes, at a cost of exactly m + 1 rotated
copies, and no odd m ever does.  So each ring offers one usable separation per
even m rather than one in total, and the candidate pool widens accordingly.

For each hub, each rational ring about it whose spindle rotation lies in the
field, and each even m, look for two ring points at the chord across m steps
and ask the exact question: is the graph uncolourable once the hub is forbidden
from sharing a colour with EITHER of them?  That is the disjunction, and it is
strictly weaker than the forced pair every earlier scan hunted.

Also reported: whether either half holds on its own, because a disjunction that
only holds when one of its halves already does is the ordinary spindle wearing
a disguise.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, width_two_spindle_separation, forced_same
from pysat.solvers import Solver

which = sys.argv[1] if len(sys.argv) > 1 else "Sa"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
MS = [2, 4, 6, 8]
t0 = time.time()
P = build_Sa(K) if which == "Sa" else build_G(K, as_graph=False)
n = len(P)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
E = sorted(set((min(a, b), max(a, b)) for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
print(f"{which}: {n} pts, {len(E)} edges, k={k}, m in {MS}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

sq_all = {}
for i in range(n):
    d = rows - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for j in np.nonzero(rat)[0]:
        j = int(j)
        if j != i:
            sq_all[(i, j)] = Fr(int(sq[j, 0]), D2)

cands = []
for u in range(n):
    ring = defaultdict(list)
    for j in range(n):
        if j != u and (u, j) in sq_all:
            ring[sq_all[(u, j)]].append(j)
    for D, mem in ring.items():
        if D == 1 or len(mem) < 2 or not closable_distance(D):
            continue
        wants = {}
        for m in MS:
            v = width_two_spindle_separation(D, m)
            if v and v > 0:
                wants.setdefault(v, m)
        if not wants:
            continue
        for a in range(len(mem) - 1):
            for b in range(a + 1, len(mem)):
                w, wp = mem[a], mem[b]
                s = sq_all.get((w, wp))
                if s in wants and (u, w) not in Eset and (u, wp) not in Eset:
                    cands.append((u, D, wants[s], w, wp))
print(f"{len(cands)} (hub, ring, m, pair) candidates  [{time.time()-t0:.0f}s]",
      flush=True)
bym = defaultdict(int)
for u, D, m, w, wp in cands:
    bym[(m, D)] += 1
print("   top:", sorted(bym.items(), key=lambda t: -t[1])[:8], flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in E:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sel0 = 1 + n * k
for idx, (u, D, m, w, wp) in enumerate(cands):
    s = sel0 + idx
    for c in range(k):
        cls.append([-s, -(1 + u * k + c), -(1 + w * k + c)])
        cls.append([-s, -(1 + u * k + c), -(1 + wp * k + c)])
sv = Solver(name="cd15", bootstrap_with=cls)
print(f"solver built, {len(cls)} clauses  [{time.time()-t0:.0f}s]", flush=True)

hits = []
for idx, (u, D, m, w, wp) in enumerate(cands):
    if not sv.solve(assumptions=[sel0 + idx]):
        half = (forced_same(sv, u, w, k), forced_same(sv, u, wp, k))
        hits.append((u, D, m, w, wp, half))
        tag = "GENUINE DISJUNCTION" if not any(half) else "one half alone"
        print(f"*** hub {u}, ring D={D}, m={m}, pair ({w},{wp}): {tag} "
              f"halves={half} ***  [{time.time()-t0:.0f}s]", flush=True)
    if idx % 5000 == 4999:
        print(f"   {idx+1}/{len(cands)}, {len(hits)} so far"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(hits)} of {len(cands)} carry it  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
