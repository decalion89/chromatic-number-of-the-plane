"""How much of G does the obstruction actually need?

G with 181 pairs forbidden does not 5-colour, and those pairs are irreducible.
But the statement is about all 1581 vertices, and there is no reason the whole
graph is load-bearing.  Removing vertices only ever makes colouring EASIER, so
any subgraph that still fails to colour is a smaller witness to the same fact
-- and a small explicit witness is what makes a result checkable by someone
else.

Two stages.  First keep only what the forbidden pairs touch plus its
neighbourhood, since a vertex far from every forbidden pair can hardly be
carrying the argument.  Then peel greedily: try each remaining vertex, drop it
if the rest still fails to colour, and repeat until a whole sweep drops
nothing.  What survives is vertex-irreducible.

The forbidden pairs are carried along: a pair whose endpoint is removed is
removed with it, so the survivor is a genuine subgraph statement and not a
bookkeeping artefact.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
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
FORB = [p for D in KEEP for p in byd[D]]
print(f"G: {n} pts, {len(E)} edges, {len(FORB)} forbidden pairs"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def colours(keep):
    kept = set(keep)
    ren = {v: i for i, v in enumerate(sorted(kept))}
    m = len(ren)
    cls = [[1 + v * k + c for c in range(k)] for v in range(m)]
    cnt = 0
    for a, b in E:
        if a in kept and b in kept:
            for c in range(k):
                cls.append([-(1 + ren[a] * k + c), -(1 + ren[b] * k + c)])
    for a, b in FORB:
        if a in kept and b in kept:
            cnt += 1
            for c in range(k):
                cls.append([-(1 + ren[a] * k + c), -(1 + ren[b] * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok, m, cnt


touched = set()
for a, b in FORB:
    touched.add(a)
    touched.add(b)
adj = defaultdict(set)
for a, b in E:
    adj[a].add(b)
    adj[b].add(a)
ring1 = set(touched)
for v in list(touched):
    ring1 |= adj[v]
for label, keep in (("all of G", set(range(n))),
                    ("touched + neighbours", ring1),
                    ("touched only", touched)):
    ok, m, cnt = colours(keep)
    print(f"  {label}: {m} vertices, {cnt} forbidden pairs kept -> "
          f"{'colours' if ok else 'DOES NOT COLOUR'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)

# Start the peel from the SMALLEST set that already fails, not from the
# neighbourhood ring: the 71 touched vertices suffice, and peeling 505 costs
# seven times as many solves per sweep for the same destination.
for start in (touched, ring1, set(range(n))):
    ok, _, _ = colours(start)
    if not ok:
        cur = set(start)
        print(f"   peeling from {len(cur)} vertices", flush=True)
        break
sweep = 0
while True:
    sweep += 1
    dropped = 0
    for v in sorted(cur):
        if v not in cur:
            continue
        trial = cur - {v}
        ok, m, cnt = colours(trial)
        if not ok:
            cur = trial
            dropped += 1
    ok, m, cnt = colours(cur)
    print(f"  sweep {sweep}: dropped {dropped}, {m} vertices left, {cnt} "
          f"forbidden pairs, colours {ok}  [{time.time()-t0:.0f}s]",
          flush=True)
    if dropped == 0:
        break
print(f"\nvertex-irreducible witness: {len(cur)} vertices of G's {n}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
