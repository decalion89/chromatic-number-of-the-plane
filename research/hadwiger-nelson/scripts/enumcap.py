"""Enumerate capped sets fast, and measure the two-distance chromatic number
of each.

Putting both conditions into one satisfiability problem was correct and far
too slow: each round ruled out a single candidate, and twenty rounds took four
hundred seconds.  The counters said why it was also unnecessary.  Over those
twenty rounds the number of sampled 5-colourings never moved, meaning every
candidate the solver proposed was genuinely capped -- capped sets are not
scarce inside the ball, they are everywhere.  What failed, every time, was the
chromatic condition.

So stop asking the solver for it.  Enumerate capped sets on their own, which
makes each round cheap -- colouring constraints plus a blocking clause per set
already seen -- and measure chi of the {1, d} graph of each afterwards.
Thousands of samples instead of twenty, over exactly the space where the
answer would have to live.

A capped set whose two-distance graph needs five colours forces a
monochromatic pair at distance d in every 5-colouring, which is the object
this whole search is for.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pairsat import pairs_at
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
M = int(sys.argv[2]) if len(sys.argv) > 2 else 12
DSQ = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(3)
HOPS = int(sys.argv[4]) if len(sys.argv) > 4 else 2
NCOL = int(sys.argv[5]) if len(sys.argv) > 5 else 200
ROUNDS = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
MINE = int(sys.argv[7]) if len(sys.argv) > 7 else 0
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E1 = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
D2 = b.D * b.D
s2 = pairs_at(b, r, int(DSQ * D2))
adjs = defaultdict(set)
for a, c in E1:
    adjs[a].add(c)
    adjs[c].add(a)
seed = pickle.load(open(SC + "core13.pkl", "rb"))[0]
cur = set(seed)
for _ in range(HOPS):
    nxt = set(cur)
    for v in cur:
        nxt |= adjs[v]
    cur = nxt
scope = sorted(cur)
ins = set(scope)
two = [(a, c) for a, c in list(E1) + list(s2) if a in ins and c in ins]
print(f"scope {len(scope)} points, {len(two)} two-distance pairs; hunting "
      f"capped {M}-sets whose {{1, sqrt{DSQ}}} graph needs {k} colours"
      f"  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E1:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve()
rng = random.Random(1000000007)


def a_colouring():
    col_solver.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                           for w in range(n * k)])
    col_solver.solve()
    mo = col_solver.get_model()
    return [next(c for c in range(k) if mo[v * k + c] > 0) for v in range(n)]


def witness_cap(T):
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


def two_chrom(T):
    loc = {v: i for i, v in enumerate(T)}
    ts = set(T)
    for kk in range(2, k + 2):
        c2 = [[1 + i * kk + c for c in range(kk)] for i in range(len(T))]
        for a, c in two:
            if a in ts and c in ts:
                for col in range(kk):
                    c2.append([-(1 + loc[a] * kk + col),
                               -(1 + loc[c] * kk + col)])
        s = Solver(name="cd15", bootstrap_with=c2)
        ok = s.solve()
        s.delete()
        if ok:
            return kk
    return k + 2


fives = [a_colouring() for _ in range(NCOL)]
print(f"seeded with {len(fives)} colourings  [{time.time()-t0:.0f}s]",
      flush=True)
pool = IDPool(start_from=n + 1)
base = list(CardEnc.equals(lits=[v + 1 for v in scope], bound=M, vpool=pool,
                           encoding=EncType.seqcounter))
# Blind enumeration drifts to sparse sets, which is the wrong half of the
# space: the chromatic number needs edges.  Demanding a floor on the
# two-distance edges inside the chosen set steers it, at the cost of one
# indicator per pair.
if MINE:
    eind = []
    for a, c in two:
        e = pool.id(("e", a, c))
        base.append([-e, a + 1])
        base.append([-e, c + 1])
        eind.append(e)
    base += list(CardEnc.atleast(lits=eind, bound=MINE, vpool=pool,
                                 encoding=EncType.seqcounter))
for si, c in enumerate(fives):
    byc = defaultdict(list)
    for v in scope:
        byc[c[v]].append(v)
    miss = [pool.id(("z", si, q)) for q in range(k)]
    base.append(miss)
    for q in range(k):
        for v in byc[q]:
            base.append([-miss[q], -(v + 1)])
solver = Solver(name="cd15", bootstrap_with=base)
best = (0, None)
seen = capped = 0
hist = defaultdict(int)
for rnd in range(ROUNDS):
    if not solver.solve():
        print(f"\nno more candidates after {seen} sets  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        break
    mod = solver.get_model()
    T = [v for v in scope if mod[v] > 0]
    solver.add_clause([-(v + 1) for v in T])       # never propose it again
    seen += 1
    w = witness_cap(T)
    if w is not None:
        solver.add_clause([pool.id(("z", len(fives), q)) for q in range(k)])
        continue
    capped += 1
    ch = two_chrom(T)
    hist[ch] += 1
    if ch > best[0]:
        best = (ch, list(T))
        print(f"   best so far: two-distance chi {ch} on a capped {M}-set "
              f"({capped} capped of {seen} proposed)  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    if ch >= k:
        pal = ring_palette_bound(cls, n * k, T, k)
        print(f"\n*** CAPPED AND {k}-CHROMATIC: {M} points, palette {pal}, "
              f"two-distance chi {ch}.  In EVERY {k}-colouring of G two of "
              f"them at distance {float(DSQ)**.5:.5f} share a colour -- "
              f"A FORCED PAIR.  [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((T, [P[v] for v in T], DSQ),
                    open(SC + f"WIN_{M}_{DSQ}.pkl".replace("/", "_"), "wb"))
        break
    if rnd % 200 == 0 and rnd:
        print(f"   {seen} proposed, {capped} capped, chi histogram "
              f"{dict(sorted(hist.items()))}  [{time.time()-t0:.0f}s]",
              flush=True)
print(f"\n{capped} capped sets of {seen} proposed; two-distance chi "
      f"histogram {dict(sorted(hist.items()))}  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
