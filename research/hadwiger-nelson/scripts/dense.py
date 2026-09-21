"""Search for a set that is capped AND dense in the second distance.

Two objects are in hand and they do not overlap.  A thirteen-point set whose
{1, sqrt3} graph needs five colours -- remarkably small, against the five
hundred a unit-distance graph needs for the same -- but whose palette in every
carrier tried is 5.  And a sixty-point capped set, palette 4, whose {1, sqrt3}
graph needs only 3.

Either alone proves nothing.  Together they would: a capped set splits into at
most four classes, each independent for distance 1, and if its {1, sqrt3}
graph needs five colours those four classes cannot all avoid sqrt3 either.
Some class then holds two points a distance sqrt3 apart, in the same colour,
in every 5-colouring.

So ask for both at once.  The decision procedure already handles the cap; the
density goes in as a cardinality constraint over the pairs at the second
distance, counted by an indicator per pair.  What comes back is capped by
construction and dense by construction, and its chromatic number for the two
distances is then one exact question.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
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
M = int(sys.argv[2]) if len(sys.argv) > 2 else 13
TARGET = int(sys.argv[3]) if len(sys.argv) > 3 else 4
MINE = int(sys.argv[4]) if len(sys.argv) > 4 else 20
DSQ = Fr(sys.argv[5]) if len(sys.argv) > 5 else Fr(3)
SEEDN = int(sys.argv[6]) if len(sys.argv) > 6 else 1200
HOPS = int(sys.argv[7]) if len(sys.argv) > 7 else 0
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = build_G(K, as_graph=False)
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
print(f"G: {n} points, {len(E1)} unit edges, {len(s2)} at distance^2 {DSQ}; "
      f"asking for {M} points, palette <= {TARGET}, at least {MINE} "
      f"two-distance edges inside  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E1:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve()
rng = random.Random(15485863)


def a_colouring():
    col_solver.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                           for w in range(n * k)])
    col_solver.solve()
    mo = col_solver.get_model()
    return [next(c for c in range(k) if mo[v * k + c] > 0) for v in range(n)]


def witness(T):
    pool = IDPool(start_from=n * k + 1)
    ind = [pool.id(("u", c)) for c in range(k)]
    ext = list(cls)
    for c in range(k):
        ext.append([-ind[c]] + [1 + v * k + c for v in T])
    ext += list(CardEnc.atleast(lits=ind, bound=TARGET + 1, vpool=pool,
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


def chrom(W, kk):
    loc = {v: i for i, v in enumerate(W)}
    ws = set(W)
    c2 = [[1 + i * kk + c for c in range(kk)] for i in range(len(W))]
    for a, c in two:
        if a in ws and c in ws:
            for col in range(kk):
                c2.append([-(1 + loc[a] * kk + col), -(1 + loc[c] * kk + col)])
    s = Solver(name="cd15", bootstrap_with=c2)
    ok = s.solve()
    s.delete()
    return ok


if HOPS:
    W13, _w = pickle.load(open(SC + "core13.pkl", "rb"))
    adjs = defaultdict(set)
    for a, c in E1:
        adjs[a].add(c)
        adjs[c].add(a)
    cur = set(W13)
    for _ in range(HOPS):
        nxt = set(cur)
        for v in cur:
            nxt |= adjs[v]
        cur = nxt
    scope = sorted(cur)
else:
    scope = list(range(n))
inscope = set(scope)
two = [(a, c) for a, c in two if a in inscope and c in inscope]
print(f"scope: {len(scope)} points, {len(two)} two-distance pairs inside",
      flush=True)
samples = [a_colouring() for _ in range(SEEDN)]
print(f"seeded with {len(samples)} colourings  [{time.time()-t0:.0f}s]",
      flush=True)
for rnd in range(40000):
    pool = IDPool(start_from=n + 1)
    y = [v + 1 for v in scope]
    f = list(CardEnc.equals(lits=y, bound=M, vpool=pool,
                            encoding=EncType.seqcounter))
    eind = []
    for a, c in two:
        e = pool.id(("e", a, c))
        eind.append(e)
        f.append([-e, a + 1])
        f.append([-e, c + 1])
    f += list(CardEnc.atleast(lits=eind, bound=MINE, vpool=pool,
                              encoding=EncType.seqcounter))
    for si, c in enumerate(samples):
        byc = defaultdict(list)
        for v in scope:
            byc[c[v]].append(v)
        miss = [pool.id(("z", si, q)) for q in range(k)]
        f += list(CardEnc.atleast(lits=miss, bound=k - TARGET, vpool=pool,
                                  encoding=EncType.seqcounter))
        for q in range(k):
            for v in byc[q]:
                f.append([-miss[q], -(v + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    ok = s.solve()
    mod = s.get_model() if ok else None
    s.delete()
    if not ok:
        print(f"\nUNSATISFIABLE after {len(samples)} colourings: no {M} "
              f"points of G are both capped at {TARGET} and carry {MINE} "
              f"two-distance edges.  An absence proof."
              f"  [{time.time()-t0:.0f}s]", flush=True)
        break
    T = [v for v in scope if mod[v] > 0]
    w = witness(T)
    if w is None:
        pal = ring_palette_bound(cls, n * k, T, k)
        ok4 = chrom(T, 4)
        ein = sum(1 for a, c in two if a in set(T) and c in set(T))
        print(f"\n*** CAPPED AND DENSE: {M} points, palette {pal}, {ein} "
              f"two-distance edges, 4-colourable for the two distances: "
              f"{ok4}  [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((T, [P[v] for v in T], DSQ),
                    open(SC + f"dense_{M}_{MINE}.pkl", "wb"))
        if not ok4:
            print(f"    *** AND NOT 4-COLOURABLE: in EVERY 5-colouring of G "
                  f"two of these points at distance {float(DSQ)**.5:.5f} "
                  f"share a colour -- A FORCED PAIR", flush=True)
            break
        MINE = ein + 2
        print(f"    raising the density demand to {MINE}", flush=True)
        continue
    samples.append(w)
    if rnd % 200 == 0 and rnd:
        print(f"   round {rnd}: {len(samples)} colourings"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
