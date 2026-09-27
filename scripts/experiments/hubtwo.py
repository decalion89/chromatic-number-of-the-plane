"""A spindle for disjunctions of width two.

The ordinary spindle needs a NAMED forced pair, and every disjunction found so
far -- 161 pairs on the witness, three antipodal pairs on Sa -- refuses it,
because each pair carries its own centre of rotation and one rotation cannot
serve them all.

Make them share the centre.  Let u be a hub, let W be a set of points all at
the same distance from u, and suppose that in every proper k-colouring u is
monochromatic with SOME point of W.  Now rho, the rotation about u joining
that ring's distance, fixes u, so every rotated copy rho^t(H) states the same
disjunction about the same u, over rho^t(W).  The colour class of u must
therefore contain one point from each rho^t(W), pairwise non-adjacent -- an
independent transversal inside a single circle.

On that circle write points as exponents of rho.  Two are at unit distance
when their exponents differ by one, by the definition of rho.  For W two
points m apart the transversal is a walk that must switch branch at every
step, and it exists for every m except m = +-2:

    W = {0, -2}:  a0 in {0,-2}, a1 in {1,-1}, a2 in {2,0}
    a0 = 0  -> both of a1 are one away.  dead.
    a0 = -2 -> a1 = 1 -> both of a2 are one away.  dead.

So three copies suffice, and the hypothesis is a distance: rho^2 moves a ring
point by 2cos(theta/2), hence |w - w'|^2 = 4 - 1/D for a ring of squared
radius D.  That is one exact test per (hub, ring, pair).

Calibrated on Sa at four colours, where de Grey's own weak property lives.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
t0 = time.time()
P = build_Sa(K)
n = len(P)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
E = sorted(set((min(a, b), max(a, b)) for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
print(f"Sa: {n} pts, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]", flush=True)

# exact squared distance between every ordered pair, as a Fraction when rational
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

# candidates: hub u, ring D about u, pair w,w' on it with |w-w'|^2 = 4 - 1/D
cands = []
for u in range(n):
    ring = defaultdict(list)
    for j in range(n):
        if j != u and (u, j) in sq_all:
            ring[sq_all[(u, j)]].append(j)
    for D, mem in ring.items():
        if D == 1 or len(mem) < 2 or not closable_distance(D):
            continue
        want = 4 - Fr(1, 1) / D
        if want <= 0:
            continue
        for a in range(len(mem) - 1):
            for b in range(a + 1, len(mem)):
                w, wp = mem[a], mem[b]
                if sq_all.get((w, wp)) == want and (u, w) not in Eset \
                        and (u, wp) not in Eset:
                    cands.append((u, D, w, wp))
print(f"{len(cands)} (hub, ring, pair) candidates at separation^2 = 4 - 1/D"
      f"  [{time.time()-t0:.0f}s]", flush=True)
byD = defaultdict(int)
for u, D, w, wp in cands:
    byD[D] += 1
print("   by ring:", dict(sorted(byD.items(), key=lambda t: -t[1])[:10]), flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in E:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sel0 = 1 + n * k
for idx, (u, D, w, wp) in enumerate(cands):
    s = sel0 + idx
    for c in range(k):
        cls.append([-s, -(1 + u * k + c), -(1 + w * k + c)])
        cls.append([-s, -(1 + u * k + c), -(1 + wp * k + c)])
sv = Solver(name="cd15", bootstrap_with=cls)
print(f"solver built, {len(cls)} clauses  [{time.time()-t0:.0f}s]", flush=True)

hits = []
for idx, (u, D, w, wp) in enumerate(cands):
    if not sv.solve(assumptions=[sel0 + idx]):
        hits.append((u, D, w, wp))
        print(f"*** WIDTH-TWO HUB DISJUNCTION: hub {u}, ring D={D}, "
              f"pair ({w},{wp}) ***  [{time.time()-t0:.0f}s]", flush=True)
    if idx % 2000 == 1999:
        print(f"   {idx+1}/{len(cands)}, {len(hits)} so far"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(hits)} of {len(cands)} carry it  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
