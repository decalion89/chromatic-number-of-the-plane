"""Is there ANY tight set at five colours, ring or not?

Every search here has been over rings, and for a good reason: a bite can only
turn a set onto itself if the turn stabilises it, and about a centre that
means a ring.  But that restriction answers a question about the MACHINE, not
about the plane.  Before asking whether a tight set can be exploited it is
worth asking whether one exists at all, and for that the set is free.

A set of m points is generically able to show all k colours.  One that cannot
-- palette at most k-1 with no centre adjacent to it to make that free -- is
the first sign of rigidity at five colours, whatever shape it has.  At four
colours they are everywhere: de Grey's is a hexagon held to two.

Searched by hill climbing.  Start from a connected seed, evaluate the palette,
and try swapping one point at a time for a nearby one, keeping whatever does
not raise it.  The palette bound is monotone under adding graph vertices but
NOT under changing the set, so this is a real search rather than a sweep --
and every evaluation is exact, so a reported tight set is proved tight.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import pickle
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

SEED = sys.argv[1] if len(sys.argv) > 1 else "G"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
M = int(sys.argv[3]) if len(sys.argv) > 3 else 7
BUDGET = int(sys.argv[4]) if len(sys.argv) > 4 else 6000
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
if SEED.endswith(".pkl"):
    P = pickle.load(open(SC + SEED, "rb"))
else:
    P = {"Sa": build_Sa, "Y": build_Y,
         "G": lambda f: build_G(f, as_graph=False)}[SEED](K)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
print(f"{SEED}: {n} points, {len(E)} edges, k={k}, sets of {M}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
cache = {}


def pal(S):
    key = tuple(sorted(S))
    if key not in cache:
        cache[key] = ring_palette_bound(cls, n * k, list(key), k)
    return cache[key]


def near(S):
    out = set()
    for v in S:
        out |= adj[v]
        for u in adj[v]:
            out |= adj[u]
    return sorted(out - set(S))


rng = random.Random(15486047)
hot = sorted(range(n), key=lambda v: -len(adj[v]))
best_overall = (k + 1, None)
evals = 0
restarts = 0
while evals < BUDGET:
    restarts += 1
    v = hot[rng.randrange(min(200, n))]
    S = {v}
    while len(S) < M:
        cand = near(S) or list(range(n))
        S.add(cand[rng.randrange(len(cand))])
    cur = pal(S)
    evals += 1
    stall = 0
    while stall < 40 and evals < BUDGET:
        out = rng.choice(sorted(S))
        cands = near(S)
        if not cands:
            break
        inn = cands[rng.randrange(len(cands))]
        T = (S - {out}) | {inn}
        if len(T) != M:
            stall += 1
            continue
        p = pal(T)
        evals += 1
        if p <= cur:
            if p < cur:
                stall = 0
            S, cur = T, p
        else:
            stall += 1
    if cur < best_overall[0]:
        best_overall = (cur, sorted(S))
        print(f"*** palette {cur} (generic is {k}) on {M} points: "
              f"{sorted(S)}  [{time.time()-t0:.0f}s]", flush=True)
        if cur <= k - 1:
            with open(SC + f"tight_{SEED}_{k}_{M}.pkl", "wb") as f:
                pickle.dump([P[i] for i in sorted(S)], f)
    if restarts % 20 == 0:
        print(f"   {restarts} restarts, {evals} evaluations, best "
              f"{best_overall[0]}  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{SEED} at k={k}: best palette on {M} points is "
      f"{best_overall[0]} after {evals} exact evaluations"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
