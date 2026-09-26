"""Many capped sets, and the two-distance chromatic number of each.

Both halves of the reduction have now been enumerated from one side.  Three
hundred and thirty-seven distinct 5-chromatic {1, sqrt3} cores live inside G,
some of them only nine points, and every one has palette 5 -- none is capped.

So try the other side.  Only ONE capped set has ever been examined here: the
seven points the decision procedure returned first, grown greedily to sixty.
Its {1, sqrt3} graph needs three colours, which is nowhere near five, but
there is no reason to think that set is typical.  Different random seeds give
the procedure different starting points, and each grows into a different
capped set.

Every subset of a capped set is capped, so a capped set whose two-distance
graph needs five colours would contain a capped core immediately -- and the
whole argument would close.  This generates capped sets and measures that
number for each.
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
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
DSQ = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(3)
TRIES = int(sys.argv[3]) if len(sys.argv) > 3 else 40
M0 = int(sys.argv[4]) if len(sys.argv) > 4 else 7
t0 = time.time()
SC = ("/tmp/hn/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E1 = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(set)
for a, c in E1:
    adj[a].add(c)
    adj[c].add(a)
dm, D2 = b.dim, b.D * b.D
tval = int(DSQ * D2)
at_d = defaultdict(set)
for u in range(n):
    d = r - r[u]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    for w in np.nonzero(ok & (sq[:, 0] == tval))[0]:
        w = int(w)
        if w != u:
            at_d[u].add(w)
print(f"G: {n} points, distance^2 {DSQ} has "
      f"{sum(len(v) for v in at_d.values())//2} pairs"
      f"  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E1:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve()


def two_chrom(W, kk):
    loc = {v: i for i, v in enumerate(W)}
    ws = set(W)
    c2 = [[1 + i * kk + c for c in range(kk)] for i in range(len(W))]
    for u in W:
        for v in (adj[u] | at_d[u]) & ws:
            if u < v:
                for col in range(kk):
                    c2.append([-(1 + loc[u] * kk + col),
                               -(1 + loc[v] * kk + col)])
    s = Solver(name="cd15", bootstrap_with=c2)
    ok = s.solve()
    s.delete()
    return ok


best = (0, None)
for trial in range(TRIES):
    rng = random.Random(7919 * (trial + 1) + 13)

    def a_colouring():
        col_solver.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                               for w in range(n * k)])
        col_solver.solve()
        mo = col_solver.get_model()
        return [next(c for c in range(k) if mo[v * k + c] > 0)
                for v in range(n)]

    def witness(T):
        pool = IDPool(start_from=n * k + 1)
        ind = [pool.id(("u", c)) for c in range(k)]
        ext = list(cls)
        for c in range(k):
            ext.append([-ind[c]] + [1 + v * k + c for v in T])
        ext += list(CardEnc.atleast(lits=ind, bound=k, vpool=pool,
                                    encoding=EncType.seqcounter))
        s = Solver(name="cd15", bootstrap_with=ext)
        ok = s.solve()
        out = None
        if ok:
            mo = s.get_model()
            out = [next(c for c in range(k) if mo[v * k + c] > 0)
                   for v in range(n)]
        s.delete()
        return out

    samples = [a_colouring() for _ in range(700)]
    found = None
    for rnd in range(60):
        pool = IDPool(start_from=n + 1)
        f = list(CardEnc.equals(lits=list(range(1, n + 1)), bound=M0,
                                vpool=pool, encoding=EncType.seqcounter))
        for si, c in enumerate(samples):
            byc = defaultdict(list)
            for v in range(n):
                byc[c[v]].append(v)
            miss = [pool.id(("z", si, q)) for q in range(k)]
            f.append(miss)
            for q in range(k):
                for v in byc[q]:
                    f.append([-miss[q], -(v + 1)])
        s = Solver(name="cd15", bootstrap_with=f)
        ok = s.solve()
        mod = s.get_model() if ok else None
        s.delete()
        if not ok:
            break
        T = [v for v in range(n) if mod[v] > 0]
        w = witness(T)
        if w is None:
            found = T
            break
        samples.append(w)
    if found is None:
        continue
    cur = list(found)
    pool_c = set()
    for v in cur:
        pool_c |= adj[v] | at_d[v]
    pool_c -= set(cur)
    grew = True
    while grew:
        grew = False
        for v in sorted(pool_c, key=lambda v: -(len((adj[v] | at_d[v])
                                                    & set(cur)))):
            if ring_palette_bound(cls, n * k, cur + [v], k) <= k - 1:
                cur.append(v)
                pool_c.discard(v)
                pool_c |= (adj[v] | at_d[v]) - set(cur)
                grew = True
                break
    ch = 2
    while ch < 8 and not two_chrom(cur, ch):
        ch += 1
    if ch > best[0]:
        best = (ch, list(cur))
    print(f"   trial {trial}: capped set of {len(cur)} points, "
          f"two-distance chi = {ch} (need 5)  [{time.time()-t0:.0f}s]",
          flush=True)
    if ch >= 5:
        print(f"\n*** CAPPED AND 5-CHROMATIC FOR {{1, sqrt{DSQ}}}: "
              f"{len(cur)} points.  In EVERY {k}-colouring of G two of them "
              f"at distance {float(DSQ)**.5:.5f} share a colour.", flush=True)
        pickle.dump((cur, [P[v] for v in cur], DSQ),
                    open(SC + f"win_{DSQ}.pkl".replace("/", "_"), "wb"))
        break
print(f"\nbest two-distance chi over the capped sets tried: {best[0]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
