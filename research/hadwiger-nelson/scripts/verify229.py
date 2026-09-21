"""Verify the 229-pair obstruction independently before it is claimed.

Greedy minimisation reports that G at five colours stops colouring once 229
pairs over five closable distance classes -- 15/16 (1 pair), 16 (12), 17/2
(24), 9 (48), 7 (144) -- are forbidden together, and that every one of the
five classes is needed.  That is the first positive statement about G at five
colours here, so it gets checked rather than reported.

Three checks.  A second and third solver on the same formula, since a single
solver's UNSAT is a single solver's.  Per-class irreducibility re-run from
scratch rather than trusted from the greedy pass, which could have dropped a
class whose necessity only appeared later.  And the pairs re-derived from the
exact geometry, so that "distance 15/16" means what it says and the obstruction
is about the plane rather than about an indexing slip.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
KEEP = [Fr(15, 16), Fr(16), Fr(17, 2), Fr(9), Fr(7)]
P = build_G(K, as_graph=False)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
n = len(P)
E = sorted(set((min(a, b), max(a, b))
               for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
byd = defaultdict(list)
for i in range(n - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in Eset:
            v = int(sq[off, 0])
            if v:
                byd[Fr(v, D2)].append((i, j))

# --- the geometry, re-derived exactly --------------------------------------
print("the five classes, distances re-derived from the exact coordinates:",
      flush=True)
allpairs = []
for D in KEEP:
    prs = byd[D]
    bad = sum(1 for i, j in prs if P[i].dist2(P[j]) != D)
    clos = closable_distance(D)
    print(f"   D={str(D):6s} {len(prs):4d} pairs, {bad} with the wrong exact "
          f"distance, closable {clos}", flush=True)
    assert bad == 0 and clos
    allpairs += prs
print(f"   {len(allpairs)} pairs in all; G has {len(E)} edges, so the "
      f"forbidden set is {100*len(allpairs)/len(E):.1f}% of them", flush=True)


def colours(pairs, solver):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in list(E) + list(pairs):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name=solver, bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


print("\nthree solvers on the full 229-pair formula:", flush=True)
for s in ("cd15", "glucose4", "minisat22"):
    ok = colours(allpairs, s)
    print(f"   {s:10s} -> {'colours' if ok else 'DOES NOT COLOUR'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)

print("\nirreducibility, each class dropped in turn and re-tested:",
      flush=True)
for D in KEEP:
    rest = [p for DD in KEEP if DD != D for p in byd[DD]]
    ok = colours(rest, "cd15")
    print(f"   without D={str(D):6s} ({len(byd[D]):4d} pairs, "
          f"{len(rest)} left) -> "
          f"{'colours -- so it IS needed' if ok else 'still does not colour'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("\nDONE", flush=True)
