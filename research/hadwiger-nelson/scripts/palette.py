"""Sweep the palette bound, which is the criterion de Grey's chain really uses.

The pattern census settled what Sa says about its hexagon at four colours: not
"some antipodal pair is monochromatic" but the far stronger "the six points
take at most TWO colours", ten surviving patterns out of a hundred and
eighty-seven.  That is a cap of k - 2, and it is why six crossing edges turn
two copies into a forced pair.  At five colours the same hexagon takes all
five and every one of the 202 patterns survives.

So the design target is a number: a ring whose palette is capped below k, and
ideally at k - 2.  Unlike the gateway it has a gradient, and unlike the forced
pair it is cheap -- one SAT call screens a ring, since a ring that can show
all k colours is capped at k and needs no further questions.

Swept over every centre and every ring of a graph, reporting the whole
distribution rather than just the minimum, so the shape of the failure is
visible.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, doubly_usable_ring, ring_palette_bound
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

which = sys.argv[1] if len(sys.argv) > 1 else "G"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
MINRING = int(sys.argv[3]) if len(sys.argv) > 3 else 6
t0 = time.time()
if which == "Sa":
    P = build_Sa(K)
elif which == "G":
    P = build_G(K, as_graph=False)
else:                                    # thickened: Sa stacked m times
    m = int(which)
    P = build_Sa(K)
    rho = rotation_joining(4, K)
    cur, seen = list(P), set(P)
    for _ in range(m):
        cur = [rho(p) for p in cur]
        for q in cur:
            if q not in seen:
                seen.add(q)
                P.append(q)
b = IntBasis.covering(P)
r = b.rows(P)
dm, d2 = b.dim, b.D * b.D
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{which}: {n} pts, {len(E)} edges, k={k}, rings of >= {MINRING}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
assert len(E) > 2 * n, "edge finder degenerated -- basis overflow"

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
nvars = n * k

rings = []
for ci in range(n):
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), d2)
        if v:
            grp[v].append(int(off))
    for D, mem in grp.items():
        if len(mem) >= MINRING and closable_distance(D):
            rings.append((ci, D, mem))
print(f"{len(rings)} (centre, ring) with >= {MINRING} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# the screen: can the ring show all k colours at once?  one call each.
ind = [nvars + 1 + c for c in range(k)]
enc = CardEnc.atleast(lits=list(ind), bound=k, top_id=nvars + k,
                      encoding=EncType.seqcounter)
sel0 = max([nvars + k] + [abs(x) for cl in enc.clauses for x in cl]) + 1
extra = list(enc.clauses)
for idx, (ci, D, mem) in enumerate(rings):
    s = sel0 + idx
    for c in range(k):
        extra.append([-s, -ind[c]] + [1 + v * k + c for v in mem])
sv = Solver(name="cd15", bootstrap_with=cls + extra)
assert sv.solve(), "the graph itself does not colour"
capped = []
for idx, (ci, D, mem) in enumerate(rings):
    if not sv.solve(assumptions=[sel0 + idx]):
        capped.append(idx)
        print(f"*** CAPPED: centre {ci}, ring D={D}, {len(mem)} points, "
              f"cannot show all {k} colours ***  [{time.time()-t0:.0f}s]",
              flush=True)
    if idx % 2000 == 1999:
        print(f"   {idx+1}/{len(rings)}, {len(capped)} capped"
              f"  [{time.time()-t0:.0f}s]", flush=True)
sv.delete()
print(f"\n{len(capped)} of {len(rings)} rings are capped below {k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for idx in capped[:20]:
    ci, D, mem = rings[idx]
    bound = ring_palette_bound(cls, nvars, mem, k)
    print(f"   centre {ci}, D={D}, {len(mem)} pts: palette {bound}, "
          f"{'DOUBLY USABLE' if doubly_usable_ring(D) else 'spindle only'}",
          flush=True)
print("DONE", flush=True)
