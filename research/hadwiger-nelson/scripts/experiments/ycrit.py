"""Is Y critical for its forcing?

The ball of radius 1.5 about the segment holds 773 of Y's 791 vertices and the
pair comes apart; the ball of radius 2.0 is all 791 and it does not.  So the
eighteen outermost vertices matter collectively.  The sharper question is
whether they matter one at a time: if Y minus ANY single vertex separates the
pair, then Y is vertex-critical for the property, and there is no smaller
forcer inside it at all.

Most queries are satisfiable and quick -- a vertex whose removal frees the
pair answers in milliseconds.  A vertex whose removal leaves it forced costs a
full unsatisfiability proof, and is the interesting case: it would mean Y has
slack.
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


def still_forced(drop):
    keep = [i for i in range(len(Y)) if i != drop]
    idx = {v: i for i, v in enumerate(keep)}
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(keep))]
    for a, b in E:
        if a in idx and b in idx:
            for c in range(4):
                cls.append([-(1 + idx[a] * 4 + c), -(1 + idx[b] * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            return "not 4-colourable"
        sv.conf_budget(60000)
        r = sv.solve_limited(assumptions=[1 + idx[ia] * 4,
                                          -(1 + idx[ib] * 4)])
    return {False: "still FORCED", True: "separates", None: "budget"}[r]


slack, hard = [], []
# A full pass is thirteen hours: separating a pair in Y minus a vertex costs
# up to 35000 conflicts (that is the range Y's own unforced pairs cost), and
# a vertex that leaves the pair forced costs a full refutation.  A stratified
# sample of sixty answers the question -- is there ANY slack -- in minutes.
allv = sorted((v for v in range(len(Y)) if v not in (ia, ib)),
              key=lambda v: deg.get(v, 0))
order = allv[::max(1, len(allv) // 60)]
for n, v in enumerate(order):
    r = still_forced(v)
    if r == "still FORCED":
        slack.append(v)
        print(f"  *** vertex {v} (degree {deg.get(v,0)}) is SLACK: Y - v "
              f"still forces  [{time.time()-t0:.0f}s]", flush=True)
    elif r == "budget":
        hard.append(v)
    if n % 10 == 0:
        print(f"  ... {n}/{len(order)}, {len(slack)} slack, {len(hard)} over "
              f"budget  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(slack)} vertices are slack, {len(hard)} inconclusive; "
      f"{'no slack in the sample -- Y looks vertex-critical' if not slack and not hard else 'Y has slack'}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
