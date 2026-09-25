"""The other direction: can few distance-avoiding sets COVER the plane?

The whole session has pushed the lower bound -- find a graph that cannot be
coloured.  The upper bound has not moved since Isbell in 1950: seven colours
by a hexagonal tiling, and nobody has ever built a six-colouring, let alone
a five.

Density says there is room.  A colour class is a set avoiding distance 1,
the densest known has density about 0.229, and five classes covering the
plane need only 1/5 = 0.2 each.  So nothing in the counting forbids a
five-colouring; what is missing is an arrangement.

Asked here on a torus, discretised, and CONSERVATIVELY: two cells get an
edge whenever some point of one is exactly a unit from some point of the
other, so a proper colouring of the cells is a genuine colouring of the
continuum.  A satisfying assignment would be a periodic colouring of the
plane with that many colours -- and unsatisfiable says nothing, because the
conservative edges constrain more than the plane does.

The test of the discretisation is whether it can reproduce Isbell: seven
colours must come out satisfiable, or the cells are too coarse to say
anything about six.
"""
import sys, time, math
from pysat.solvers import Solver

SIDE = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
CELLS = int(sys.argv[2]) if len(sys.argv) > 2 else 30
t0 = time.time()
d = SIDE / CELLS
print(f"torus {SIDE} x {SIDE}, {CELLS}x{CELLS} cells of side {d:.4f} "
      f"(diagonal {d*math.sqrt(2):.4f})  [{time.time()-t0:.0f}s]", flush=True)
n = CELLS * CELLS


def cell_centre(i):
    return ((i % CELLS + 0.5) * d, (i // CELLS + 0.5) * d)


def torus_gap(i, j):
    """Smallest and largest distance between any two points of cells i, j,
    measured on the torus."""
    xi, yi = cell_centre(i)
    xj, yj = cell_centre(j)
    dx = abs(xi - xj)
    dx = min(dx, SIDE - dx)
    dy = abs(yi - yj)
    dy = min(dy, SIDE - dy)
    lo = math.hypot(max(0.0, dx - d), max(0.0, dy - d))
    hi = math.hypot(dx + d, dy + d)
    return lo, hi


E = []
for i in range(n):
    for j in range(i + 1, n):
        lo, hi = torus_gap(i, j)
        if lo <= 1.0 <= hi:
            E.append((i, j))
print(f"{n} cells, {len(E)} conservative edges, mean degree "
      f"{2*len(E)/n:.1f}  [{time.time()-t0:.0f}s]", flush=True)
for KC in (7, 6, 5):
    cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
    for a, b in E:
        for c in range(KC):
            cls.append([-(1 + a * KC + c), -(1 + b * KC + c)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    note = ""
    if KC == 7 and not ok:
        note = "   <- the discretisation is too coarse to say anything"
    if ok and KC < 7:
        note = "   *** a periodic colouring of the plane with " \
               f"{KC} colours ***"
    print(f"   {KC} colours: {'SAT' if ok else 'UNSAT'}{note}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
