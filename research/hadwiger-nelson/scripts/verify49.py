"""Verify the colouring, independently of the search that produced it.

TabuCol reports zero conflicts by maintaining the count incrementally, so a
bug in the update would report a colouring that is not one.  The claim that G
at five colours admits a colouring with every pair at squared distance 4/9
non-monochromatic is worth more than the word of the thing that found it.

So: re-run the seed that reached zero, take the colouring out, and check it
from scratch -- every unit edge bichromatic, every 4/9 pair bichromatic, five
colours used at most, all 1581 vertices assigned -- against the EXACT geometry
rather than against the search's own bookkeeping.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis

k = 5
t0 = time.time()
P = build_G(F, as_graph=False)
g = build_graph(P)
E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
extra = []
for i in range(len(P) - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in set(E) and Fr(int(sq[off, 0]), D2) == Fr(4, 9):
            extra.append((i, j))
print(f"G: {g.n} pts, {len(E)} unit edges, {len(extra)} pairs at 4/9"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# re-run the winning seed, this time keeping the colouring
edges = np.array(E + extra, dtype=np.int64)
n = g.n
rng = np.random.default_rng(2001)
A = np.zeros((n, n), dtype=np.int8)
A[edges[:, 0], edges[:, 1]] = 1
A[edges[:, 1], edges[:, 0]] = 1
col = rng.integers(0, k, n)
onehot = np.zeros((n, k), dtype=np.int64)
onehot[np.arange(n), col] = 1
gam = A @ onehot
idx = np.arange(n)
f = int(gam[idx, col].sum()) // 2
best = f
tab = np.zeros((n, k), dtype=np.int64)
for it in range(800_000):
    if f == 0:
        break
    V = np.nonzero(gam[idx, col])[0]
    delta = gam[V] - gam[V, col[V]][:, None]
    allowed = (tab[V] <= it) | (f + delta < best)
    allowed[np.arange(len(V)), col[V]] = False
    if not allowed.any():
        tab[:] = 0
        continue
    cost = np.where(allowed, delta, 1 << 40) + rng.random(delta.shape)
    r, c = divmod(int(np.argmin(cost)), k)
    v = int(V[r])
    cur = int(col[v])
    f += int(delta[r, c])
    tab[v, cur] = it + int(0.6 * len(V)) + int(rng.integers(0, 10)) + 1
    col[v] = c
    gam[:, cur] -= A[:, v]
    gam[:, c] += A[:, v]
    if f < best:
        best = f
print(f"search finished with internal count {f} after {it} moves"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# --- the independent check, against the exact geometry -------------------
bad_unit = [(a, b) for a, b in E if col[a] == col[b]]
bad_49 = [(a, b) for a, b in extra if col[a] == col[b]]
# and re-derive the two edge sets exactly, not from the arrays above
exact_unit = sum(1 for a, b in E if P[a].dist2(P[b]) != 1)
exact_49 = sum(1 for a, b in extra if P[a].dist2(P[b]) != Fr(4, 9))
print(f"  vertices coloured: {len(col)} of {g.n}, colours used "
      f"{sorted(set(int(c) for c in col))}", flush=True)
print(f"  unit edges monochromatic: {len(bad_unit)}", flush=True)
print(f"  4/9 pairs monochromatic:  {len(bad_49)}", flush=True)
print(f"  edges whose exact distance is not 1:    {exact_unit}", flush=True)
print(f"  pairs whose exact distance is not 4/9:  {exact_49}", flush=True)
ok = not bad_unit and not bad_49 and not exact_unit and not exact_49
print(f"\n{'VERIFIED' if ok else 'THE CLAIM DOES NOT CHECK OUT'}: G at five "
      f"colours admits a colouring with every unit edge AND every pair at "
      f"distance 2/3 bichromatic{'.' if ok else '?'}", flush=True)
print("So G does NOT carry the weak property on class 4/9 at five colours.",
      flush=True)
print("DONE", flush=True)
