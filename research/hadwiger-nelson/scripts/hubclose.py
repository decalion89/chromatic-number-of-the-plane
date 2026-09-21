"""Close G about its hub, which is what Sa is.

Sa carries the weak property on two rings; S, its 39-point seed, carries none
and has no rational ring about the origin at all.  What the closure adds is
not size but a CENTRE: Sa is the twelve-element dihedral closure of S about
the origin, and the origin is Sa's hub -- degree 60, the highest in the graph.
Symmetry and density at the same point.

G has a hub too.  G[0] is the image of that origin under build_G's rotation
and it also has degree 60, and G's rational ring structure is richest there --
60, 48, 36, 36, 24, 12, 12 points on D = 1, 5/3, 1/3, 4/3, 5/9, 4, 3.  But
G's SYMMETRY is about (-2, 0), the spindle pivot, which is not a hub.  Sa has
the two coincident and G does not.

So make them coincide: close G dihedrally about G[0] and ask the gateway
question there.  That is the direct analogue of what Sa is, one level up.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
G = build_G(F, as_graph=False)
C = G[0]
g0 = build_graph(G)
deg = [0] * len(G)
for a, b in g0.edges():
    deg[a] += 1
    deg[b] += 1
print(f"G: {len(G)} pts; G[0] degree {deg[0]}, max degree {max(deg)}, "
      f"{sum(1 for d in deg if d == max(deg))} vertices at the max"
      f"  [{time.time()-t0:.0f}s]", flush=True)

rot = _rot60(F).about(C)
seen, U = set(), []
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, C.y + (C.y - q.y))
            if q not in seen:
                seen.add(q)
                U.append(q)
print(f"closed about the hub: {len(U)} pts  [{time.time()-t0:.0f}s]",
      flush=True)

g = build_graph(U)
E = set((min(a, b), max(a, b)) for a, b in g.edges())
idx = {p: i for i, p in enumerate(U)}
by = defaultdict(list)
for i, p in enumerate(U):
    d = p.dist2(C)
    if not all(x == 0 for x in d.c[1:]) or d.c[0] == 0:
        continue
    q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
    j = idx.get(q)
    if j is not None and i < j and (i, j) not in E:
        by[d.c[0]].append((i, j))
order = sorted(by, key=lambda d: -len(by[d]))
print(f"graph: {g.n} pts, {len(E)} edges, {len(order)} rational rings about "
      f"the hub  [{time.time()-t0:.0f}s]", flush=True)
print(f"  biggest: {[(str(d), len(by[d])) for d in order[:8]]}", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
for a, b in sorted(E):
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sel0 = g.n * k + 1
for t, D in enumerate(order):
    for i, j in by[D]:
        for c in range(k):
            cls.append([-(1 + i * k + c), -(1 + j * k + c), sel0 + t])
sv = Solver(name="cd15", bootstrap_with=cls)
print(f"solving ({len(cls)} clauses)  [{time.time()-t0:.0f}s]", flush=True)
if not sv.solve():
    print("*** NOT 5-COLOURABLE -- chi >= 6 ***", flush=True)
else:
    hits = 0
    for t, D in enumerate(order):
        if not sv.solve(assumptions=[-(sel0 + t)]):
            hits += 1
            print(f"  *** WEAK PROPERTY: D={D}, {len(by[D])} antipodal pairs "
                  f"at squared distance {4*D}, spindle closable "
                  f"{closable_distance(4*D)} ***  [{time.time()-t0:.0f}s]",
                  flush=True)
        if t % 20 == 0:
            print(f"    ring {t}/{len(order)}, {hits} hits"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"{hits} of {len(order)} rings carry the weak property"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
