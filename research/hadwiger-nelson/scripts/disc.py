"""Cut by geometry before cutting by solver: the smallest disc that still refuses four.

Greedy deletion pays about forty-six seconds for every vertex it removes, since
a removal has to be certified by a refutation, so shedding three hundred
vertices is four hours.  But refusing four is monotone in the vertex set, and
the vertices of a graph sit in the PLANE -- so ask a geometric question first:

    what is the smallest disc about c whose points still refuse four?

Monotone in the radius, so binary search finds it in about a dozen refutations
instead of hundreds, and it can shed a large outer shell in one answer.  A
critical subgraph has no reason to be spatially compact, but it has no reason
not to be either, and a dozen solves is cheap enough to find out.

Tried about several centres -- the graph's centroid, its densest vertex, and a
few random ones -- since the disc that works depends entirely on where the
strain actually sits.
"""
import sys, time, json, math, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 4
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247.json"
CENTRES = int(sys.argv[2]) if len(sys.argv) > 2 else 6
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
X = lambda v, c: 1 + v * K + c

def refuses(keep):
    ks = set(keep)
    if len(ks) < 10: return False
    cnf = [[X(v, c) for c in range(K)] for v in keep]
    for x, y in E:
        if x in ks and y in ks:
            for c in range(K):
                cnf.append([-X(x, c), -X(y, c)])
    tri = None
    for u in keep:
        for v in sorted(adj[u] & ks):
            w = adj[u] & adj[v] & ks
            if w: tri = (u, v, min(w)); break
        if tri: break
    if tri:
        for i, v in enumerate(tri):
            cnf.append([X(v, i)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    r = s.solve(); s.delete()
    return not r

print(f"{NAME} n={n} edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)
cx0 = sum(hx) / n; cy0 = sum(hy) / n
dense = max(range(n), key=lambda v: len(adj[v]))
random.seed(0)
centres = [("centroid", cx0, cy0), ("densest", hx[dense], hy[dense])]
for k in range(CENTRES - 2):
    v = random.randrange(n)
    centres.append((f"vertex {v}", hx[v], hy[v]))
best = (n, None)
for name, cx, cy in centres:
    order = sorted(range(n), key=lambda v: (hx[v]-cx)**2 + (hy[v]-cy)**2)
    lo, hi = 10, n                     # hi is known to refuse four
    if not refuses(order[:hi]):
        print(f"  {name}: whole graph does not refuse four?", flush=True); break
    while lo < hi:
        mid = (lo + hi) // 2
        if refuses(order[:mid]): hi = mid
        else: lo = mid + 1
    keep = order[:hi]
    m = sum(1 for x, y in E if x in set(keep) and y in set(keep))
    r = math.sqrt((hx[keep[-1]]-cx)**2 + (hy[keep[-1]]-cy)**2)
    print(f"  {name:<12s}: smallest disc keeps {hi} vertices, {m} edges, "
          f"radius {r:.3f}   [{time.time()-t0:.0f}s]", flush=True)
    if hi < best[0]:
        best = (hi, keep)
        json.dump({"source": NAME, "centre": name, "n": hi,
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in g.vertices[v].x.c],
                               [[t.numerator, t.denominator] for t in g.vertices[v].y.c]]
                              for v in keep]},
                  open(f"{ROOT}/data/disc_{hi}.json", "w"))
        print(f"    written data/disc_{hi}.json", flush=True)
print(f"\n  best disc: {best[0]} vertices   [{time.time()-t0:.0f}s]", flush=True)
