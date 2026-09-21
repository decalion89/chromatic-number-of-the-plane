"""Grow the capped set along the two-distance structure, not along the index.

The forced colour class inside a capped set is independent in the unit-distance
graph, and it escapes the argument by also avoiding the chosen distance d.  So
the obstacle is the largest set that is independent for BOTH distances at
once -- the independence number of the two-distance graph {1, d} restricted to
the capped set -- and the argument closes exactly when that number drops below
|T| / palette(T).

Which is where the cap earns its keep.  With no cap the bar is |T|/5 and the
two-distance graph would have to be 6-chromatic.  With palette 4 the bar is
|T|/4 and 5-chromatic is enough: one whole step cheaper, and the step is the
one nobody can climb.

The first capped set was grown in index order, so its two-distance structure
is whatever happened to come along.  This grows it deliberately: among the
points that keep the palette at 4, take the one that adds the most edges at
distance 1 or d.  Same exact verification at every step, same upward
monotonicity -- what holds here holds in every supergraph.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 4
DSQ = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(3)
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
Eset = set(E)
n = len(P)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
dm, D2 = b.dim, b.D * b.D


def d2(u, v):
    d = r[u] - r[v]
    sq = (b._field_square(d[None, :dm]) + b._field_square(d[None, dm:]))[0]
    return Fr(int(sq[0]), D2) if not any(sq[1:]) else None


# every point at distance^2 DSQ from each vertex, once
at_d = defaultdict(set)
sq_all = None
import numpy as np
for u in range(n):
    d = r - r[u]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    tgt = int(DSQ * D2)
    for w in np.nonzero(ok & (sq[:, 0] == tgt))[0]:
        at_d[u].add(int(w))
print(f"distance^2 {DSQ}: {sum(len(v) for v in at_d.values())//2} pairs in G"
      f"  [{time.time()-t0:.0f}s]", flush=True)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)


def links(v, cur):
    s = set(cur)
    return len(adj[v] & s) + len(at_d[v] & s)


def alpha_escapes(T):
    """Largest set independent for both distances; compare to the pigeonhole."""
    m = len(T)
    q = -(-m // TARGET)
    pool = IDPool(start_from=m + 1)
    f = list(CardEnc.atleast(lits=list(range(1, m + 1)), bound=q, vpool=pool,
                             encoding=EncType.seqcounter))
    for a in range(m):
        for c in range(a + 1, m):
            u, v = T[a], T[c]
            if (min(u, v), max(u, v)) in Eset or v in at_d[u]:
                f.append([-(a + 1), -(c + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    ok = s.solve()
    s.delete()
    return ok, q


seed, _ = pickle.load(open(SC + f"grown_{k}_{TARGET}.pkl", "rb"))
cur = sorted(seed[:7])
print(f"seed: {cur}, palette {ring_palette_bound(cls, n*k, cur, k)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
pool_c = set()
for v in cur:
    pool_c |= adj[v] | at_d[v]
    for u in list(adj[v] | at_d[v]):
        pool_c |= adj[u] | at_d[u]
pool_c -= set(cur)
best = None
while True:
    cands = sorted(pool_c, key=lambda v: (-links(v, cur), v))
    took = None
    for v in cands:
        t = sorted(cur + [v])
        if ring_palette_bound(cls, n * k, t, k) <= TARGET:
            took = v
            cur = t
            break
    if took is None:
        break
    pool_c.discard(took)
    pool_c |= (adj[took] | at_d[took]) - set(cur)
    esc, q = alpha_escapes(cur)
    if not esc:
        print(f"*** BLOCKS: {len(cur)} points, forced class >= {q}, and no "
              f"class of that size avoids distance^2 {DSQ}.  In EVERY "
              f"5-colouring of G some pair at distance {float(DSQ)**.5:.5f} "
              f"is monochromatic -- A FORCED PAIR.  [{time.time()-t0:.0f}s]",
              flush=True)
        with open(SC + f"forcedpair_{DSQ}.pkl", "wb") as fh:
            pickle.dump((cur, [P[v] for v in cur], DSQ), fh)
        break
    if len(cur) % 20 == 0:
        dd = sum(1 for i, u in enumerate(cur) for v in cur[i + 1:]
                 if v in at_d[u])
        print(f"   {len(cur)} points, forced class {q}, {dd} pairs at "
              f"distance^2 {DSQ}  [{time.time()-t0:.0f}s]", flush=True)
    if len(cur) > 400:
        break
print(f"\nfinal: {len(cur)} points  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
