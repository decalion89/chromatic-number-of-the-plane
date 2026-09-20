"""de Grey's shape of forcing, asked directly.

A single forced pair is too much to hope for -- the necklace has none and the
neighbourhood theorem says forcing is never local.  What de Grey actually uses
is weaker and disjunctive: for a vertex u, in EVERY 4-colouring at least one
point at distance sqrt3 from u shares u's colour.  Call that property (*) at
u.  It is one SAT call:

    is there a 4-colouring with c(v) != c(u) for every v at sqrt3 from u?

UNSAT means (*) holds.  And sqrt3 is exactly the rhombus tip distance, the one
the Moser rotation spindles, which is why that distance and not another.

Run over every vertex of Sa, with the size of each sqrt3-sphere recorded --
the smaller the sphere on which (*) holds, the tighter the gadget, and the
fewer copies a pigeonhole argument would need afterwards.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.graph import build_graph
from pysat.solvers import Solver

pts = build_Sa()
g = build_graph(pts)
E = list(g.edges())
n = len(pts)
print(f"Sa: {n} points, {len(E)} edges", flush=True)

three = [[] for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        dx = pts[j].x - pts[i].x
        dy = pts[j].y - pts[i].y
        if dx * dx + dy * dy == 3:
            three[i].append(j)
sizes = sorted({len(s) for s in three})
print(f"sqrt3-sphere sizes present: {sizes}", flush=True)


def x(v, c):
    return 1 + v * 4 + c


cls = [[x(v, c) for c in range(4)] for v in range(n)]
for v in range(n):
    for c in range(4):
        for d in range(c + 1, 4):
            cls.append([-x(v, c), -x(v, d)])
for a, b in E:
    for c in range(4):
        cls.append([-x(a, c), -x(b, c)])

t0 = time.time()
solver = Solver(name="cd19", bootstrap_with=cls)
hits = []
try:
    for u in range(n):
        S = three[u]
        if not S:
            continue
        # c(u) = 0 without loss, and no v in S may take colour 0
        ass = [x(u, 0)] + [-x(v, 0) for v in S]
        if not solver.solve(assumptions=ass):
            hits.append((u, len(S)))
            print(f"  *** (*) HOLDS at vertex {u}, sqrt3-sphere of "
                  f"{len(S)}  [{time.time()-t0:.0f}s]", flush=True)
        if u % 50 == 0:
            print(f"  ... {u}/{n}, {len(hits)} so far "
                  f"[{time.time()-t0:.0f}s]", flush=True)
finally:
    solver.delete()
print(f"{len(hits)} vertices where every 4-colouring puts u's colour on its "
      f"sqrt3-sphere  [{time.time()-t0:.0f}s]", flush=True)
if hits:
    print(f"  smallest sphere with the property: "
          f"{min(h[1] for h in hits)}", flush=True)
