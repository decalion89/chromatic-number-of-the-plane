"""Five points of one unit circle, pairwise forced different: the apex route.

If five points at distance 1 from a common centre are forced to take five
different colours in every proper k-colouring, then the centre has k colours
in its neighbourhood and cannot be coloured: the cone over them needs k+1.
No spindle, no disjunction, no rotation -- one configuration and one
conclusion.

The geometry constrains it sharply.  A neighbourhood lies on a unit circle,
where two points are a unit apart only if they subtend 60 degrees, so the
graph a circle carries is a subgraph of disjoint hexagons: bipartite, with
maximum clique 2.  At most TWO of the five can be pairwise adjacent, and
every other pair must be forced different by the carrier while not being
adjacent at all.  So the question is how much non-adjacent forcing the
carrier supplies inside a neighbourhood, and a five-clique in the
forced-different relation is the whole target.

Forced-different is one solve: by colour permutation, if u and v CAN share a
colour they can share colour 0, so the test is two assumptions.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, itertools
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "G"
KC = int(sys.argv[2]) if len(sys.argv) > 2 else 5
LIM = int(sys.argv[3]) if len(sys.argv) > 3 else 400
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, c in E:
    for col in range(KC):
        cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve(), f"{CAR} is not {KC}-colourable"
print(f"{CAR}: {n} pts, {len(E)} edges, k = {KC}; a {KC}-clique of "
      f"forced-different points on one circle would need {KC+1} colours"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def forced_diff(u, v):
    return not s.solve(assumptions=[1 + u * KC, 1 + v * KC])


order = sorted(range(n), key=lambda v: -len(adj[v]))
bestclique = 0
extra_total = 0
checked = 0
for ci, v in enumerate(order[:LIM]):
    N = sorted(adj[v])
    if len(N) < KC:
        continue
    checked += 1
    # the forced-different graph on N(v): adjacency is automatic, the rest
    # has to be earned
    fd = defaultdict(set)
    extra = 0
    for u, w in itertools.combinations(N, 2):
        if w in adj[u]:
            fd[u].add(w)
            fd[w].add(u)
        elif forced_diff(u, w):
            fd[u].add(w)
            fd[w].add(u)
            extra += 1
    extra_total += extra
    # largest clique in fd, by simple growth from each vertex
    best = 0
    for u in N:
        cl2 = [u]
        cand = set(fd[u])
        while cand:
            w = max(cand, key=lambda x: len(fd[x] & cand))
            cl2.append(w)
            cand &= fd[w]
        best = max(best, len(cl2))
    if best > bestclique:
        bestclique = best
        print(f"   v{v} (deg {len(N)}): largest forced-different clique "
              f"{best}, {extra} non-adjacent forced pairs"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if best >= KC:
        print(f"   *** v{v} has {KC} pairwise forced-different neighbours: "
              f"its cone needs {KC+1} colours ***", flush=True)
        break
    if ci % 50 == 49:
        print(f"   {ci+1} centres, best clique {bestclique}, "
              f"{extra_total} non-adjacent forced pairs total"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{checked} neighbourhoods examined; largest forced-different clique "
      f"{bestclique} (need {KC}); {extra_total} non-adjacent forced pairs "
      f"found in total  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
