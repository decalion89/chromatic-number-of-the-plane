"""Close the bitten graph under the symmetry the seed was closed under.

De Grey's chain is one level of a recipe, not a one-off.  A 39-point seed is
closed under the twelve-element dihedral group about the origin to give Sa,
397 points, which carries a RING-shaped weak property at four colours; the
bite turns that into Y, 791 points, which forces a named pair; the spindle
turns that into G.

Running the recipe one level up means treating the level-one output as the new
seed and closing it under the same group.  Sa is already closed, so the orbit
of Y is Sa together with the twelve images of the bitten half -- about five
thousand points, which is exactly the scale the measured factor of ten
predicts for five colours, and which the hub closure of G (11047 points) came
in above.

Tested with the gateway at five: every centre, every ring whose antipodal
pairs have a spindle, forbidden all at once in one SAT call.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Sb, build_Y
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, doubly_usable_ring
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
rot60 = _rot60(K)
Y = build_Y(K)
seen, U = set(), []
for p in Y:
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                U.append(q)
            q = rot60(q)
print(f"dihedral closure of Y: {len(Y)} -> {len(U)} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)

b = IntBasis.covering(U)
r = b.rows(U)
dm, d2 = b.dim, b.D * b.D
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(U)
print(f"{n} pts, {len(E)} edges, {len(E)/n:.2f} per vertex"
      f"  [{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
sv0 = Solver(name="cd15", bootstrap_with=cls)
base_ok = sv0.solve()
print(f"{k}-colourable: {base_ok}  [{time.time()-t0:.0f}s]", flush=True)
if not base_ok:
    print(f"*** NOT {k}-COLOURABLE -- chi >= {k+1} ***", flush=True)
    sys.exit()
sv0.delete()

key = {tuple(r[i]): i for i in range(n)}
cands = []
for ci in range(n):
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for m in range(1, dm):
        ok &= sq[:, m] == 0
    ring = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), d2)
        if v:
            ring[v].append(int(off))
    for D, mem in ring.items():
        ms = set(mem)
        pairs = []
        for j in mem:
            anti = key.get(tuple(2 * r[ci] - r[j]))
            if anti is not None and anti in ms and anti > j:
                pairs.append((j, anti))
        if len(pairs) >= 2 and closable_distance(4 * D):
            cands.append((ci, D, pairs))
    if ci % 1000 == 999:
        print(f"   centres {ci+1}/{n}, {len(cands)} candidates"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"{len(cands)} (centre, ring) candidates  [{time.time()-t0:.0f}s]",
      flush=True)
sizes = defaultdict(int)
for ci, D, pairs in cands:
    sizes[len(pairs)] += 1
print("   antipodal pairs per candidate:", dict(sorted(sizes.items())),
      flush=True)

extra, sel0 = [], 1 + n * k
for idx, (ci, D, pairs) in enumerate(cands):
    s = sel0 + idx
    for a, c in pairs:
        for col in range(k):
            extra.append([-s, -(1 + a * k + col), -(1 + c * k + col)])
sv = Solver(name="cd15", bootstrap_with=cls + extra)
assert sv.solve()
hits = 0
for idx, (ci, D, pairs) in enumerate(cands):
    if not sv.solve(assumptions=[sel0 + idx]):
        hits += 1
        kind = "DOUBLY USABLE" if doubly_usable_ring(D) else "spindle only"
        print(f"*** GATEWAY: centre {ci}, ring D={D}, {len(pairs)} antipodal "
              f"pairs, {kind} ***  [{time.time()-t0:.0f}s]", flush=True)
    if idx % 2000 == 1999:
        print(f"   tested {idx+1}/{len(cands)}, {hits} gateways"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{hits} of {len(cands)} carry the gateway at {k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
