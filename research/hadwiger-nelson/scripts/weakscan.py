"""Does G at five colours have the weak property on ANY ring, at ANY centre?

The mechanism turned out to be two-stage, and only the second stage is the
strong hypothesis everything here has been testing.

  Stage one, weak:   in every 4-colouring of Sa, at least ONE of the three
                     antipodal pairs on the D=4 ring is monochromatic.  Sa
                     alone, no bite.  ONE SAT call: forbid all three and the
                     graph stops colouring.
  Stage two, strong: the bite rho_4 sharpens "one of three" to the named pair
                     (-2,0),(2,0) -- which Y forces, measured, UNSAT in 210s.
  Then spindle at squared distance 16 and get G.

So the gateway is the weak property, and it costs a single call.  Every scan
in this file has been hunting stage two directly, which is the expensive end.

Scan it properly: every vertex of G as centre, every rational ring about it
carrying at least two antipodal pairs, at five colours.  The ring is worth
having only if 4D is closable -- otherwise there is no spindle to follow with
-- so those go first, but the rest are cheap enough to take too.

One solver throughout.  Each ring gets a selector variable s and clauses
(-x_ic or -x_jc or s) for every antipodal pair and colour, so assuming -s
switches that ring's constraints on and the solver keeps everything it has
learned between rings.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
P = build_G(F, as_graph=False)
n = len(P)
g = build_graph(P)
E = set((min(a, b), max(a, b)) for a, b in g.edges())
idx = {p: i for i, p in enumerate(P)}
print(f"G: {n} pts, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D

# --- every (centre, ring) with at least two antipodal pairs ----------------
cands = []
for ci in range(n):
    C = P[ci]
    d = rows - rows[ci]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    by = defaultdict(list)
    for off in np.nonzero(rat)[0]:
        v = int(sq[off, 0])
        if not v:
            continue
        p = P[off]
        q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
        j = idx.get(q)
        if j is not None and off < j and (off, j) not in E:
            by[Fr(v, D2)].append((int(off), j))
    for D, pr in by.items():
        if len(pr) >= 2:
            cands.append((closable_distance(4 * D), len(pr), D, ci, pr))
    if ci % 400 == 0:
        print(f"  centre {ci}/{n}, {len(cands)} candidate rings"
              f"  [{time.time()-t0:.0f}s]", flush=True)

cands.sort(key=lambda r: (not r[0], -r[1]))
print(f"\n{len(cands)} (centre, ring) candidates; "
      f"{sum(1 for c in cands if c[0])} have a closable spindle"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# --- one solver, one selector per candidate --------------------------------
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in sorted(E):
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sel0 = n * k + 1
for t, (_, _, _, _, pr) in enumerate(cands):
    s = sel0 + t
    for i, j in pr:
        for c in range(k):
            cls.append([-(1 + i * k + c), -(1 + j * k + c), s])
sv = Solver(name="cd15", bootstrap_with=cls)
assert sv.solve(), "G is 5-colourable"
print(f"solver built: {len(cls)} clauses, {len(cands)} selectors"
      f"  [{time.time()-t0:.0f}s]", flush=True)

hits = []
for t, (spin, m, D, ci, pr) in enumerate(cands):
    if not sv.solve(assumptions=[-(sel0 + t)]):
        hits.append((D, ci, m, spin))
        print(f"  *** WEAK PROPERTY HOLDS: centre {ci}, ring D={D}, "
              f"{m} antipodal pairs at squared distance {4*D}, "
              f"spindle closable {spin} ***  [{time.time()-t0:.0f}s]",
              flush=True)
    if t and t % 500 == 0:
        print(f"  {t}/{len(cands)}, {len(hits)} hits"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nDONE: {len(hits)} of {len(cands)} rings carry the weak property"
      f"  [{time.time()-t0:.0f}s]", flush=True)
