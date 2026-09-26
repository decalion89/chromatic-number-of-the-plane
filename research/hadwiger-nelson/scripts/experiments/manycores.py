"""Enumerate 5-chromatic two-distance cores and ask each whether it is capped.

The reduction is settled: a set whose {1, d} graph needs five colours, and
whose palette inside the carrier is at most four, forces a monochromatic pair
at distance d in every 5-colouring.  The four classes a capped set is split
into are independent for distance 1 by construction, and if they were also
independent for d the partition would 4-colour a graph that needs five.

Both halves exist separately.  A thirteen-point core with chi = 5 for
{1, sqrt3} -- strikingly small next to the five hundred a single distance
needs -- and a sixty-point capped set.  The first has palette 5, the second
has chi 3.

But there is not one core.  The unsatisfiable core a solver returns depends on
the order it sees the clauses, so shuffling the graph yields different ones,
and only one of them has to be capped.  This spins that wheel: permute, pull a
core, measure its palette, keep going.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
DSQ = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(3)
TRIES = int(sys.argv[3]) if len(sys.argv) > 3 else 300
CARRIER = sys.argv[4] if len(sys.argv) > 4 else "G"
t0 = time.time()
SC = ("/tmp/hn/")
P = (build_G(K, as_graph=False) if CARRIER == "G"
     else pickle.load(open(SC + CARRIER, "rb")))
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E1 = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
dm, D2 = b.dim, b.D * b.D
tval = int(DSQ * D2)
s2 = set()
for u in range(n):
    d = r - r[u]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    for w in np.nonzero(ok & (sq[:, 0] == tval))[0]:
        w = int(w)
        if w != u:
            s2.add((min(u, w), max(u, w)))
two = sorted(set(E1) | s2)
print(f"{CARRIER}: {n} points, {len(E1)} unit edges, {len(s2)} at "
      f"distance^2 {DSQ}, {len(two)} together  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E1:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
kk, SEL = 4, n * 4
rng = random.Random(1000003)
seen = set()
best = (k + 1, None)
for t in range(TRIES):
    order = list(range(n))
    rng.shuffle(order)
    pos = {v: i for i, v in enumerate(order)}
    ed = [(pos[a], pos[c]) for a, c in two]
    rng.shuffle(ed)
    keep = set(range(n))
    for rnd in range(10):
        cl = []
        for v in keep:
            cl.append([-(SEL + 1 + v)] + [1 + v * kk + c for c in range(kk)])
        for a, c in ed:
            if a in keep and c in keep:
                for col in range(kk):
                    cl.append([-(1 + a * kk + col), -(1 + c * kk + col)])
        s = Solver(name="cd15", bootstrap_with=cl)
        ok = s.solve(assumptions=[SEL + 1 + v for v in sorted(keep)])
        core = None if ok else s.get_core()
        s.delete()
        if ok:
            keep = None
            break
        nxt = {abs(l) - SEL - 1 for l in core} if core else set(keep)
        if len(nxt) >= len(keep):
            break
        keep = nxt
    if keep is None:
        continue
    W = tuple(sorted(order[i] for i in keep))
    if W in seen:
        continue
    seen.add(W)
    pal = ring_palette_bound(cls, n * k, list(W), k)
    if pal < best[0]:
        best = (pal, W)
    if pal <= k - 1:
        print(f"\n*** CAPPED 5-CHROMATIC CORE: {len(W)} points, palette "
              f"{pal} of {k}, for distance^2 {DSQ}.  In EVERY {k}-colouring "
              f"two of them at distance {float(DSQ)**.5:.5f} share a colour."
              f"  [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((list(W), [P[v] for v in W], DSQ),
                    open(SC + f"cappedcore_{DSQ}.pkl".replace("/", "_"), "wb"))
        break
    if t % 25 == 0:
        print(f"   {t} tries, {len(seen)} distinct cores, sizes seen "
              f"{sorted({len(x) for x in seen})[:6]}, best palette "
              f"{best[0]}  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(seen)} distinct cores, best palette {best[0]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
