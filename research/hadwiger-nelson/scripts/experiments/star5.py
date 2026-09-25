"""The same lemma one colour up -- which is the shot at six.

Extracted from Y at four colours: in every 4-colouring, at least one of the
12 points at distance sqrt3 from the origin shares the origin's colour.  That
is de Grey's forcing, and sqrt3 is the distance the Moser rotation spindles,
since |1 - rho|^2 = 1/3 means the angle arccos(5/6) and |v - rho v| =
sqrt3 . (1/sqrt3) = 1.

The analogue at five colours would do the same job again:

    is there a vertex u of G such that EVERY 5-colouring puts u's colour
    somewhere on u's sqrt3-sphere?

One SAT call per vertex -- fix c(u) = 0 and forbid colour 0 on the whole
sphere; UNSAT is the property.  G is 5-chromatic so 5-colourings exist and the
question is not vacuous.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from pysat.solvers import Solver

t0 = time.time()
g = build_G()
pts = list(g.vertices)
E = list(g.edges())
n = len(pts)
print(f"G: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

# The exact comparison over a degree-16 field costs 2.5 million field
# operations here, so the spheres are found by a float grid first and only
# the survivors are checked exactly -- the same prefilter the graph builder
# itself uses.
import math
R3 = math.sqrt(3.0)
cells = {}
for i, q in enumerate(pts):
    cells.setdefault((math.floor(q.fx / R3), math.floor(q.fy / R3)),
                     []).append(i)
three = [[] for _ in range(n)]
for i, q in enumerate(pts):
    cx, cy = math.floor(q.fx / R3), math.floor(q.fy / R3)
    for a in (-2, -1, 0, 1, 2):
        for b in (-2, -1, 0, 1, 2):
            for j in cells.get((cx + a, cy + b), ()):
                if j == i:
                    continue
                d2 = (q.fx - pts[j].fx) ** 2 + (q.fy - pts[j].fy) ** 2
                if abs(d2 - 3.0) > 1e-7:
                    continue
                dx = pts[j].x - pts[i].x
                dy = pts[j].y - pts[i].y
                if dx * dx + dy * dy == 3:
                    three[i].append(j)
sz = sorted({len(s) for s in three})
print(f"sqrt3-sphere sizes: {sz}  [{time.time()-t0:.0f}s]", flush=True)


def x(v, c):
    return 1 + v * 5 + c


cls = [[x(v, c) for c in range(5)] for v in range(n)]
for a, b in E:
    for c in range(5):
        cls.append([-x(a, c), -x(b, c)])

solver = Solver(name="cd19", bootstrap_with=cls)
hits = []
try:
    order = sorted(range(n), key=lambda v: -len(three[v]))
    for k, u in enumerate(order):
        S = three[u]
        if not S:
            continue
        if not solver.solve(assumptions=[x(u, 0)] + [-x(v, 0) for v in S]):
            hits.append((u, len(S)))
            print(f"  *** (*) HOLDS AT FIVE COLOURS, vertex {u}, sphere "
                  f"{len(S)}  [{time.time()-t0:.0f}s]", flush=True)
        if k % 25 == 0:
            print(f"  ... {k}/{n} (sphere {len(S)}), {len(hits)} so far "
                  f"[{time.time()-t0:.0f}s]", flush=True)
finally:
    solver.delete()
print(f"{len(hits)} vertices with the property at five colours  "
      f"[{time.time()-t0:.0f}s]", flush=True)
