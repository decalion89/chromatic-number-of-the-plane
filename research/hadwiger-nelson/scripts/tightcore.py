"""What makes Sa tight at four?  Extract the core by UNSAT, not by deletion.

Sa refuses every circular clique below 4, which is what tightness means, and
today's synthesis says caps live exactly where a carrier is tight.  So the
structure responsible for Sa's tightness is the template for the thing that
is missing at five, and it is worth seeing rather than guessing at.

Delta-debugging paid a full refusal proof per trial and crawled from 397
points to 137 in an hour.  Activation literals do it properly: give every
vertex a literal a(v) with the clause a(v) -> v takes a position, solve
under assumptions {a(v) : v in S}, and the solver returns an UNSAT core that
is a smaller S.  Iterating converges in a handful of rounds, and what comes
out is a subgraph that still refuses the ratio -- still tight at four, and
as small as the proof needs.
"""
import sys, time, pickle
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "Sa"
R = Fr(*map(int, (sys.argv[2] if len(sys.argv) > 2 else "35/9").split("/")))
t0 = time.time()
P = {"Sa": build_Sa, "Y": build_Y,
     "G": lambda k: build_G(k, as_graph=False)}[CAR](K)
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), R.numerator, R.denominator
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)


def X(v, j):
    return 1 + v * p + j


def A(v):
    return 1 + n * p + v


cls = [[-A(v)] + [X(v, j) for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-X(a, j), -X(c, (j + d) % p)])
s = Solver(name="cd15", bootstrap_with=cls)
print(f"{CAR}: {n} pts, {len(E)} edges; refusing K({p}/{q}) = {float(R):.4f}?"
      f"  [{time.time()-t0:.0f}s]", flush=True)
keep = set(range(n))
for rnd in range(30):
    ok = s.solve(assumptions=[A(v) for v in sorted(keep)])
    if ok:
        print(f"   round {rnd}: {len(keep)} points MAP -- stopping",
              flush=True)
        break
    core = {lit - 1 - n * p for lit in s.get_core()}
    print(f"   round {rnd}: refuses, core {len(core)} of {len(keep)}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if len(core) >= len(keep):
        break
    keep = core
    pickle.dump([P[v] for v in sorted(keep)],
                open(SC + f"tight_{CAR}_{p}_{q}.pkl", "wb"))
ks = sorted(keep)
m = sum(1 for v in ks for u in adj[v] if u in keep) // 2
degs = sorted(len(adj[v] & keep) for v in ks)
import math
xy = [(float(P[v].x), float(P[v].y)) for v in ks]
diam = max(math.hypot(a[0] - c[0], a[1] - c[1]) for a in xy for c in xy)
print(f"\ntight core: {len(ks)} points, {m} edges, degrees "
      f"{degs[0]}..{degs[-1]}, diameter {diam:.3f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
pickle.dump([P[v] for v in ks], open(SC + f"tight_{CAR}_{p}_{q}.pkl", "wb"))
print("DONE", flush=True)
