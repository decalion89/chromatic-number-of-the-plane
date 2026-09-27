"""Chromatic number of lattice two-distance graphs, by SAT on growing patches.

L = Z[w] (w = e^{2 pi i/3}) scaled so that the edge set is {|z|^2 in NORMS}.
{1, 2/sqrt3}: scale 1/sqrt3, unit = norm 3, 2/sqrt3 = norm 4.
"""
import sys, time, itertools
from pysat.solvers import Solver
t0 = time.time()
def norm(x, y): return x * x - x * y + y * y        # |x + y w|^2
def test(norms, R, K):
    pts = [(x, y) for x in range(-R, R + 1) for y in range(-R, R + 1) if norm(x, y) <= R * R]
    idx = {p: i for i, p in enumerate(pts)}
    S = [(a, b) for a in range(-4, 5) for b in range(-4, 5) if norm(a, b) in norms]
    E = set()
    for (x, y), i in idx.items():
        for a, b in S:
            j = idx.get((x + a, y + b))
            if j is not None and j > i: E.add((i, j))
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(len(pts))]
    for i, j in E:
        for c in range(K): cl.append([-X(i, c), -X(j, c)])
    cl.append([X(0, 0)])
    s = Solver(name="cd19", bootstrap_with=cl); r = s.solve()
    return len(pts), len(E), r
for name, norms in [("{1,2/sqrt3} (norms 3,4)", {3, 4}), ("{1,2} (norms 1,4)", {1, 4}),
                    ("{1,sqrt3} (norms 1,3)", {1, 3}), ("{1,sqrt7/..}(3,7)", {3, 7}),
                    ("{1,2/sqrt3,sqrt7/sqrt3}(3,4,7)", {3, 4, 7})]:
    for K in (4, 5, 6):
        for R in (4, 6, 9, 13):
            n, e, r = test(norms, R, K)
            if not r:
                print(f"  {name}: NOT {K}-colourable at radius {R} (n={n}, e={e})  [{time.time()-t0:.0f}s]", flush=True)
                break
        else:
            print(f"  {name}: {K}-colourable up to radius 13  [{time.time()-t0:.0f}s]", flush=True)
            break
