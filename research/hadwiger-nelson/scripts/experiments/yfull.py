"""Is Y vertex-critical for its forcing?  All 789, no budgets.

Every vertex sampled so far is essential: deleting it lets (2,0) and (-2,0)
take different colours.  Without a budget each query costs a few seconds --
the budget was the thing making this look expensive -- so the full statement
is an hour's work rather than a day's, and it is worth having exactly.

If every vertex is essential then Y is a MINIMAL four-colour forcer for its
pair, with nothing to trim, which is the cleanest thing that can be said about
de Grey's middle floor and the sharpest reason to expect the five-colour
analogue to be large.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Y
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
Y = build_Y()
g = build_graph(Y)
E = list(g.edges())
ia = next(i for i, p in enumerate(Y) if p == Point(F.rational(2), F.zero()))
ib = next(i for i, p in enumerate(Y) if p == Point(F.rational(-2), F.zero()))
deg = {}
for a, b in E:
    deg[a] = deg.get(a, 0) + 1
    deg[b] = deg.get(b, 0) + 1
print(f"Y: {len(Y)} vertices, {len(E)} edges, pair {ia},{ib}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

slack, iso = [], []
order = [v for v in range(len(Y)) if v not in (ia, ib)]
for n, v in enumerate(order):
    if deg.get(v, 0) == 0:
        iso.append(v)
        continue
    keep = [i for i in range(len(Y)) if i != v]
    idx = {w: i for i, w in enumerate(keep)}
    cls = [[1 + w * 4 + c for c in range(4)] for w in range(len(keep))]
    for a, b in E:
        if a in idx and b in idx:
            for c in range(4):
                cls.append([-(1 + idx[a] * 4 + c), -(1 + idx[b] * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        sv.solve()
        r = sv.solve(assumptions=[1 + idx[ia] * 4, -(1 + idx[ib] * 4)])
    if not r:
        slack.append(v)
        print(f"  *** vertex {v} (degree {deg.get(v,0)}) is SLACK: Y - v "
              f"still forces  [{time.time()-t0:.0f}s]", flush=True)
    if n % 50 == 0:
        print(f"  ... {n}/{len(order)}, {len(slack)} slack so far  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(slack)} slack vertices, {len(iso)} isolated; "
      f"{'Y IS VERTEX-CRITICAL for its forcing' if not slack else 'Y has slack'}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
if slack:
    print(f"  the slack ones: {slack}", flush=True)
