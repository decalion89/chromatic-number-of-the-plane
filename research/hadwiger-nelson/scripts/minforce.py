"""Is the forcing local, and how small a graph does it need?

Y has 791 vertices and forces (2,0),(-2,0).  Nothing says it needs all of
them.  A smaller forcer would be a better building block than Y -- the whole
architecture is "forcer, then spindle", so shrinking the forcer shrinks
everything above it, and it is the object to port to another field.

Two questions, cheapest first.  Does the forcing survive inside a ball around
the segment joining the pair?  And does it survive under greedy deletion of
whole layers?  A ball that still forces is the answer; a ball that does not
says the forcing is genuinely spread and Y is doing all of the work.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Y
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
Y = build_Y()
gy = build_graph(Y)
E = list(gy.edges())
zs = [(float(p.x), float(p.y)) for p in Y]
A = Point(F.rational(2), F.zero())
B = Point(F.rational(-2), F.zero())
ia = next(i for i, p in enumerate(Y) if p == A)
ib = next(i for i, p in enumerate(Y) if p == B)
print(f"Y: {len(Y)} vertices, {len(E)} edges; the pair is {ia}, {ib}  "
      f"[{time.time()-t0:.0f}s]", flush=True)


def forced(keep):
    idx = {v: i for i, v in enumerate(sorted(keep))}
    n = len(idx)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        if a in idx and b in idx:
            for c in range(4):
                cls.append([-(1 + idx[a] * 4 + c), -(1 + idx[b] * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            return "not 4-colourable", n
        sv.conf_budget(4000000)
        r = sv.solve_limited(assumptions=[1 + idx[ia] * 4,
                                          -(1 + idx[ib] * 4)])
    return {False: "FORCED", True: "separable", None: "budget"}[r], n


# distance from the segment joining the pair, which runs along the x-axis
def dist(i):
    x, y = zs[i]
    cx = min(2.0, max(-2.0, x))
    return ((x - cx) ** 2 + y * y) ** 0.5


for r in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 10.0):
    keep = {i for i in range(len(Y)) if dist(i) <= r + 1e-9}
    if ia not in keep or ib not in keep:
        continue
    verdict, n = forced(keep)
    print(f"  ball of radius {r}: {n} vertices -- {verdict}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if verdict == "FORCED":
        print(f"  *** the forcing is local: {n} vertices suffice ***",
              flush=True)
        import pickle
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/minforce.pkl", "wb") as fh:
            pickle.dump(sorted(keep), fh)
        break
