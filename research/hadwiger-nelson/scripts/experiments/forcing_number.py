"""The one invariant everything reduces to: the smallest rainbow-forcing set.

A core of p is a set T with min |c(N(p) union T)| = k -- a set that uses every
colour in every colouring. Strip the pivot away and what is left is a single
graph invariant:

    rho(W, k) = the least size of a set using all k colours in every
                k-colouring.

Every bound in this package is a statement about it. In a k-vertex-critical
graph rho = n: for any u, colour W - u with k-1 and give u the kth, and a set
omitting u misses that colour, so every vertex is needed. That IS the
inertness of de Grey's G. In a uniquely k-colourable graph rho = k, one
representative per class -- the triangular lattice. And a core at a pivot is
a forcing set containing N(p), so the smallest blockable core, three, needs
rho at most 3 + the pivot's circle.

It is computed the same way cores are: ask whether some colouring leaves the
first colour unused on S, and if so add a vertex carrying it. One solve per
point.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import pickle
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Sa, build_Y
from hn.field import Field
from hn.forced import ColourRelations, min_colours_on
from hn.geometry import Point
from hn.graph import build_graph

SC = "/tmp/hn"


def forcing_set(rel, limit=400):
    """Smallest set using all k colours, by counterexamples."""
    S = []
    for _ in range(limit):
        rel.calls += 1
        if not rel._s.solve(assumptions=[-rel._x(v, 0) for v in S]):
            return S, True
        m = set(rel._s.get_model())
        cand = [v for v in range(rel.graph.n)
                if v not in S and rel._x(v, 0) in m]
        if not cand:
            return S, False
        S.append(min(cand, key=lambda v: -len(rel.graph.adj[v])))
    return S, False


def shrink(rel, S):
    cur, changed = list(S), True
    while changed and len(cur) > 1:
        changed = False
        for t in list(cur):
            trial = [x for x in cur if x != t]
            if trial and min_colours_on(rel, trial) >= rel.k:
                cur, changed = trial, True
                break
    return cur


# the triangular lattice, uniquely 3-colourable
f = Field((3,))
r3, half = f.sqrt(3), f.rational(Fr(1, 2))
steps = [Point(f.rational(1), f.zero()), Point(half, r3 * half),
         Point(-half, r3 * half)]
lat = {Point(f.zero(), f.zero())}
for _ in range(3):
    new = set()
    for q in lat:
        for s in steps:
            new.add(Point(q.x + s.x, q.y + s.y))
            new.add(Point(q.x - s.x, q.y - s.y))
    lat |= new
lattice = build_graph(sorted(lat, key=lambda q: float(q.x * q.x + q.y * q.y)))

mult, (i, j) = pickle.load(open(f"{SC}/topvecs.pkl", "rb"))[0]
G = build_G()
dx = G.vertices[j].x - G.vertices[i].x
dy = G.vertices[j].y - G.vertices[i].y
union = build_graph(list(dict.fromkeys(
    list(G.vertices) + [Point(v.x + dx, v.y + dy) for v in G.vertices])))

t0 = time.time()
for name, g, k in (("triangular lattice", lattice, 3),
                   ("Sa", build_Sa(), 4), ("Y", build_Y(), 4),
                   ("G", G, 5), ("G union (G+t)", union, 5)):
    g = g if hasattr(g, "vertices") else build_graph(g)
    rel = ColourRelations(g, k)
    if not rel.colourable:
        rel.close()
        continue
    S, ok = forcing_set(rel)
    if ok:
        S = shrink(rel, S)
    print(f"{name:>20} n={g.n:5d} k={k}: rho = "
          f"{len(S) if ok else '>' + str(len(S))}  "
          f"(k = {k}, n = {g.n})  [{time.time()-t0:.0f}s]", flush=True)
    rel.close()
