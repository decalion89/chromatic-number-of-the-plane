"""Is any unit-distance graph here UNIQUELY 4-colourable?  rho = k would say so.

rho = k exactly when a graph is uniquely k-colourable, since k vertices can
force all k colours only if the classes themselves are pinned.  The triangular
lattice does it at three: rho = 3, its only 3-colouring is the coset one.

And that is not a curiosity.  Two copies of the lattice at the Moser angle are
4-CHROMATIC -- which this package first met as a vacuity bug, the union
silently ceasing to be 3-colourable.  Read forwards it is a mechanism: destroy
the unique colouring and the chromatic number rises.

So the question worth asking is whether the same trick has a fourth floor. A
uniquely 4-colourable unit-distance graph would give chi >= 5 by the same
move, with the obvious analogue at five.

This tests rho = 4 on every 4-chromatic graph in the package, and checks the
lattice mechanism explicitly.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from fractions import Fraction
from hn import degrey
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from hn.mixed import three_hexagon_gadget
from rhotight import run
from pysat.solvers import Solver

FLD = Field((3, 11))
HALF = FLD.rational(Fraction(1, 2))
MOSER = Rotation(FLD.rational(Fraction(5, 6)),
                 FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))


def colourable(g, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


def tri_patch(r):
    u = Point(HALF, FLD.sqrt(3) * HALF)
    one = Point(FLD.one(), FLD.zero())
    return [one.scaled(a) + u.scaled(b)
            for a in range(-r, r + 1) for b in range(-r, r + 1)]


print("the mechanism, stated as a measurement:", flush=True)
for r in (2, 3, 4):
    base = tri_patch(r)
    g = build_graph(base)
    seen, both = set(base), list(base)
    for p in base:
        q = MOSER(p)
        if q not in seen:
            seen.add(q)
            both.append(q)
    w = build_graph(both)
    print(f"  patch radius {r}: n={g.n:4}, 3-colourable {colourable(g, 3)}; "
          f"with a Moser-rotated copy n={w.n:4}, 3-colourable "
          f"{colourable(w, 3)}", flush=True)

print("\nis any 4-chromatic graph here uniquely 4-colourable (rho = 4)?",
      flush=True)
field, pivot, gad = three_hexagon_gadget()
for name, pts in (("three-hexagon gadget", gad), ("de Grey S", degrey.build_S()),
                  ("de Grey Sa", degrey.build_Sa())):
    g = build_graph(pts)
    if not colourable(g, 4) or colourable(g, 3):
        print(f"  {name}: not 4-chromatic, skipped", flush=True)
        continue
    t = time.time()
    v = run(g, 4, 4, rounds=100000, report=10 ** 9)
    print(f"  {name}: n={g.n:4}, rho = 4 is {v}"
          + ("   *** UNIQUELY 4-COLOURABLE ***" if v is True else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
