"""Where four-colour forcing actually lives: Sa.

The necklace has none -- measured, 0 of its 8385 non-adjacent pairs -- and the
reason is structural: the neighbourhood of any vertex of a unit-distance graph
is bipartite, so no colour is ever forced locally.  Forcing has to be
long-range, and de Grey gets it out of Sa, the 397-point hexagonal object.

So look there.  For each non-adjacent pair, ask whether some 4-colouring
separates them.  Colour symmetry lets that be one SAT call with two
assumptions: fix c(u) = 0 and demand c(v) = 1, since any separating colouring
can be permuted into that form.  UNSAT means the pair is monochromatic in
EVERY 4-colouring.

If they exist, their distances say which rotation would spindle them -- and by
3(4d^2-1) being a square in F, which fields can carry the same step.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.graph import build_graph
from pysat.solvers import Solver

pts = build_Sa()
g = build_graph(pts)
E = list(g.edges())
n = len(pts)
print(f"Sa: {n} points, {len(E)} edges", flush=True)
adj = {}
for a, b in E:
    adj.setdefault(a, set()).add(b)
    adj.setdefault(b, set()).add(a)


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
forced = []
try:
    for i in range(n):
        for j in range(i + 1, n):
            if j in adj.get(i, ()):
                continue
            if solver.solve(assumptions=[x(i, 0), x(j, 1)]):
                continue
            forced.append((i, j))
        if i % 40 == 0:
            print(f"  ... vertex {i}/{n}, {len(forced)} forced pairs "
                  f"[{time.time()-t0:.0f}s]", flush=True)
finally:
    solver.delete()
print(f"{len(forced)} pairs monochromatic in every 4-colouring  "
      f"[{time.time()-t0:.0f}s]", flush=True)

if forced:
    dist = Counter()
    for i, j in forced:
        dx = pts[j].x - pts[i].x
        dy = pts[j].y - pts[i].y
        dist[str(dx * dx + dy * dy)] += 1
    print("squared distances of the forced pairs:", flush=True)
    for d, c in dist.most_common(12):
        print(f"    {d}   x{c}", flush=True)
    import pickle
    with open("/tmp/claude-0/-home-user-darwin-50/"
              "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/saforced.pkl",
              "wb") as fh:
        pickle.dump(forced, fh)
