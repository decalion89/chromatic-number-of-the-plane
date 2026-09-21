"""One call: is Y's antipodal pair (-2,0),(2,0) forced at four colours?

The ring criterion predicted it and the geometry confirms it -- the pivot
build_G rotates about is in Y, exactly one vertex of Y sits at squared
distance 16 from it, and the composition of build_G's two rotations is
cos 31/32, sin sqrt(63)/32, which is the spindle rotation for D = 16 exactly.
So the prediction is that adding the edge between them makes Y uncolourable
at four.  That is one SAT call, and if it comes back UNSAT the template is
established rather than inferred.

The controls matter as much: the same pair at FIVE colours must be free (Y is
comfortably 5-colourable), and Sa alone -- Y's half, before the biting
rotation -- must not force it either, or the rotation is doing nothing.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
A = Point(F.rational(-2), F.zero())
B = Point(F.rational(2), F.zero())


def forced(P, name, k):
    idx = {p: i for i, p in enumerate(P)}
    if A not in idx or B not in idx:
        print(f"{name}: pair not present (A={A in idx}, B={B in idx})",
              flush=True)
        return
    i, j = idx[A], idx[B]
    g = build_graph(P)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"{name} at {k}: NOT {k}-COLOURABLE", flush=True)
        return
    free = sv.solve(assumptions=[1 + i * k, -(1 + j * k)])
    print(f"{name} at {k} colours ({g.n} pts, {len(E)} edges): "
          f"pair {i},{j} is {'FREE' if free else 'FORCED SAME'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()


Y = build_Y(F)
forced(Y, "Y", 4)
forced(Y, "Y", 5)
forced(build_Sa(F), "Sa (before the biting rotation)", 4)
print("DONE", flush=True)
