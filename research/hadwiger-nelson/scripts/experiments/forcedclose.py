"""Is there ANY pair of points the carrier forces to within q/2?

This is the open end, stated as a search rather than a hope.  The spindle
closes only when two points are forced to within half the adjacency
threshold -- at K(9/2) that means positions equal or adjacent, distance 0 or
1.  De Grey's bite produces exactly this in colours (forced-same) and, by
complete enumeration of 9^7 assignments, produces the opposite in positions:
a floor on the spread, not a ceiling.

So ask the carrier directly.  For a pair (u, v), "forced to within 1" means
no homomorphism of the whole graph puts them 2 or more apart, which is one
UNSAT per pair: assert that their circular distance is at least 2 and see if
anything survives.  If some pair is forced, a spindle about one of them
turns it into a contradiction and chi_c climbs past 9/2.  If no pair is,
that is the mechanism's death certificate in this carrier, written in
UNSATs rather than guessed at.

Pairs are taken in order of geometric distance, nearest first: two points a
long way apart have no reason to be tied, and the ones an edge could later
be spindled onto are the close ones.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, math
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "Sa"
R = Fr(*map(int, (sys.argv[2] if len(sys.argv) > 2 else "9/2").split("/")))
NPAIR = int(sys.argv[3]) if len(sys.argv) > 3 else 200
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), R.numerator, R.denominator
KMAX = q // 2                      # the spindle needs forcing to within this
xy = [(float(t.x), float(t.y)) for t in P]
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve(), f"{CAR} does not map to K({p}/{q})"
print(f"{CAR}: {n} pts, {len(E)} edges, K({p}/{q}); a spindle needs a pair "
      f"forced to within {KMAX}  [{time.time()-t0:.0f}s]", flush=True)
# nearest non-adjacent pairs first
cand = []
for u in range(n):
    for v in range(u + 1, n):
        if v in adj[u]:
            continue
        d = math.hypot(xy[u][0] - xy[v][0], xy[u][1] - xy[v][1])
        if d < 2.0:
            cand.append((d, u, v))
cand.sort()
print(f"{len(cand)} non-adjacent pairs closer than 2 apart; testing the "
      f"nearest {min(NPAIR, len(cand))}  [{time.time()-t0:.0f}s]", flush=True)
SEL = n * p + 1
best = None
for k, (d, u, v) in enumerate(cand[:NPAIR]):
    # assert circular distance >= KMAX+1 under a selector
    sel = SEL + k
    for j in range(p):
        allowed = [jj for jj in range(p)
                   if min((jj - j) % p, (j - jj) % p) > KMAX]
        s.add_clause([-sel, -(1 + u * p + j)] + [1 + v * p + jj
                                                 for jj in allowed])
    if not s.solve(assumptions=[sel]):
        print(f"   *** pair {u},{v} at distance {d:.4f} is FORCED to within "
              f"{KMAX} ***  [{time.time()-t0:.0f}s]", flush=True)
        best = (u, v, d)
        break
    if k % 25 == 24:
        print(f"   {k+1} pairs tested, none forced  [{time.time()-t0:.0f}s]",
              flush=True)
if best is None:
    print(f"\nno pair among the nearest {min(NPAIR,len(cand))} is forced to "
          f"within {KMAX}: the spindle has nothing to act on in {CAR}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
else:
    u, v, d = best
    print(f"\nFORCED PAIR {u},{v} at distance {d:.6f} -- spindling about "
          f"either endpoint would close it", flush=True)
    pickle.dump((u, v), open(SC + f"forced_{CAR}_{p}_{q}.pkl", "wb"))
print("DONE", flush=True)
