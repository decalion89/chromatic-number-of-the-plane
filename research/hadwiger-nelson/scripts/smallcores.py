"""Small 5-chromatic cores, and the palette of each -- the cheap side of the
same question.

Asking a solver for a capped set carrying as many two-distance edges as a
5-chromatic core needs is the expensive way round: the cardinality bound is
tight and the first call runs for minutes.  The cheap way is to take the cores
themselves, which cost seconds to extract, and measure their palettes.  A core
of nine points carries nineteen edges by construction; if one of them is also
capped, the argument closes.

Cores are pulled from the ball with the clause order shuffled, so each call
returns a different one, and only the small sizes matter -- nine to eleven
points, where the edge threshold is low enough that a capped set might reach
it.
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
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
DSQ = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(3)
TRIES = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
MAXSZ = int(sys.argv[4]) if len(sys.argv) > 4 else 12
HOPS = int(sys.argv[5]) if len(sys.argv) > 5 else 2
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
print(f"scope {len(scope)} points, {len(two)} two-distance pairs; cores up "
      f"to {MAXSZ} points  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E1:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
kk, SEL = 4, n * 4
rng = random.Random(2038074743)
seen = set()
hist = defaultdict(lambda: defaultdict(int))
best = (k + 1, None)
for t in range(TRIES):
    ed = list(two)
    rng.shuffle(ed)
    keep = set(scope)
    for rnd in range(12):
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
    if keep is None or len(keep) > MAXSZ:
        continue
    W = tuple(sorted(keep))
    if W in seen:
        continue
    seen.add(W)
    pal = ring_palette_bound(cls, n * k, list(W), k)
    hist[len(W)][pal] += 1
    if pal < best[0]:
        best = (pal, W)
    if pal <= k - 1:
        print(f"\n*** CAPPED 5-CHROMATIC CORE: {len(W)} points, palette "
              f"{pal}.  In EVERY {k}-colouring of G two of them at distance "
              f"{float(DSQ)**.5:.5f} share a colour -- A FORCED PAIR."
              f"  [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((list(W), [P[v] for v in W], DSQ),
                    open(SC + f"WINCORE_{DSQ}.pkl".replace("/", "_"), "wb"))
        break
    if t % 300 == 0 and t:
        print(f"   {t} tries, {len(seen)} distinct small cores, palettes by "
              f"size { {s: dict(v) for s, v in sorted(hist.items())} }, best "
              f"{best[0]}  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(seen)} distinct cores of at most {MAXSZ} points; palettes by "
      f"size { {s: dict(v) for s, v in sorted(hist.items())} }; best "
      f"{best[0]}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
