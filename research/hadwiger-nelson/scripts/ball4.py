"""How much graph does one forced pair cost?  Measured by growing a ball.

At four colours the ceiling is cheap to find: five_247_c is 5-chromatic and
vertex-critical, so H = G - p is 4-colourable and mu_4(H, p) = 4 for every one
of its 803 vertices.  Pairwise forcing needs the cap attained -- |N(p)| = 4 --
and exactly two vertices qualify.  At p = 315, N = {130, 461, 561, 757} holds
one 60-degree edge, so five pairs there are genuinely forced apart.

The question is what one costs.  UNSAT-core extraction on the whole 802-vertex
instance is the wrong end to start from: proving that graph 4-uncolourable is
already hard, and the selector encoding makes it harder.  Growing a ball is the
right end.  Order the vertices by graph distance from the pair, and binary
search the smallest PREFIX in which the pair is still forced apart.  Every
small prefix is a small instance and answers instantly; the search costs about
ten solves, and the answer is a number.

A few dozen vertices means a forced pair is a LOCAL object, and the failure to
find one at five is a failure of search.  Essentially all of H means it is a
symptom of being one vertex short of (k+1)-chromatic, and hunting one as a
stepping stone to exactly that is circular.
"""
import sys, time, json
from fractions import Fraction as Fr
from itertools import combinations
from collections import deque
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 4
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = [set() for _ in range(n)]
for x, y in E:
    adj[x].add(y); adj[y].add(x)
print(f"five_247_c n={n} m={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

def forced_apart(keep, u, v, budget=None):
    idx = {w: i for i, w in enumerate(keep)}
    X = lambda w, c: 1 + idx[w] * K + c
    cnf = [[X(w, c) for c in range(K)] for w in keep]
    ks = set(keep)
    for x, y in E:
        if x in ks and y in ks:
            for c in range(K):
                cnf.append([-X(x, c), -X(y, c)])
    cnf.append([X(u, 0)]); cnf.append([X(v, 0)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if budget is None:
        r = s.solve()
    else:
        s.conf_budget(budget); r = s.solve_limited()
    s.delete()
    if r is None: return None
    return not r

PIVOT = 315
nb = sorted(adj[PIVOT])
alive = [w for w in range(n) if w != PIVOT]
aset = set(alive)
print(f"pivot {PIVOT}, N = {nb}, |N| = {len(nb)} = the cap", flush=True)

def order_from(u, v):
    """vertices of H by graph distance from {u, v}, nearest first"""
    dist = {u: 0, v: 0}; q = deque([u, v]); out = [u, v]
    while q:
        w = q.popleft()
        for z in adj[w]:
            if z in aset and z not in dist:
                dist[z] = dist[w] + 1; q.append(z); out.append(z)
    out += [w for w in alive if w not in dist]      # anything unreachable
    return out, dist

report = {}
for u, v in combinations(nb, 2):
    d2 = float((g.vertices[u] - g.vertices[v]).norm2())
    edge = v in adj[u]
    full = forced_apart(alive, u, v, budget=2_000_000)
    if full is not True:
        print(f"  ({u},{v}) d^2={d2:.6f} edge={edge}  forced={full}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        continue
    order, dist = order_from(u, v)
    lo, hi = 2, len(order)                    # hi is known to work
    while lo < hi:
        mid = (lo + hi) // 2
        if forced_apart(order[:mid], u, v, budget=2_000_000) is True:
            hi = mid
        else:
            lo = mid + 1
    pref = order[:hi]
    rad = max(dist.get(w, 99) for w in pref)
    m = sum(1 for x, y in E if x in set(pref) and y in set(pref))
    print(f"  ({u},{v}) d^2={d2:.6f} edge={edge}  FORCED; smallest prefix "
          f"{hi} vertices ({m} edges), radius {rad}   [{time.time()-t0:.0f}s]",
          flush=True)
    report[f"{u},{v}"] = {"d2": d2, "edge": edge, "prefix": hi,
                          "edges": m, "radius": rad, "vertices": sorted(pref)}
json.dump({"pivot": PIVOT, "neighbourhood": nb, "prefixes": report},
          open(f"{ROOT}/data/forced_certificates.json", "w"), indent=1)
print(f"\nwritten   [{time.time()-t0:.0f}s]", flush=True)
