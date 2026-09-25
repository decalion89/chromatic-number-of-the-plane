"""What makes rho small?  Measure it exactly, wherever that is affordable.

The decision on de Grey's G is out of reach, so the useful question is the
design one: across the graphs this package can measure exactly, which
structural features put rho near its floor k and which push it to n?

Two extremes are already known. A uniquely k-colourable graph has rho = k --
the triangular lattice at three colours, whose only 3-colouring is the coset
one. A k-vertex-critical graph has rho = n, since colouring G - u with k-1 and
giving u the kth leaves u alone carrying it.

Everything else sits between, and the table is the design knowledge: it says
what a construction has to look like for rho to be small enough that a core of
three exists.

Each value is decided, not deleted: budget B gives a forcing set or is
refuted, and the smallest B that still gives one is rho.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
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


def chi(g, hi=6):
    for k in range(2, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
        for a, b in g.edges():
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return None


def tri_patch(r=4):
    u = Point(HALF, FLD.sqrt(3) * HALF)
    one = Point(FLD.one(), FLD.zero())
    return [one.scaled(a) + u.scaled(b)
            for a in range(-r, r + 1) for b in range(-r, r + 1)]


def moser():
    u = Point(HALF, FLD.sqrt(3) * HALF)
    one = Point(FLD.one(), FLD.zero())
    rh = [Point(FLD.zero(), FLD.zero()), one, u, one + u]
    rot = Rotation(FLD.rational(Fraction(5, 6)),
                   FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
    return rh + [rot(p) for p in rh[1:]]


field, pivot, gad = three_hexagon_gadget()
CASES = [("triangular patch", tri_patch()), ("Moser spindle", moser()),
         ("three-hexagon gadget", gad), ("de Grey S", degrey.build_S()),
         ("de Grey Sa", degrey.build_Sa())]

print(f"{'graph':22} {'n':>6} {'chi':>4} {'rho':>5}  reading")
print("-" * 62)
for name, pts in CASES:
    g = build_graph(pts)
    k = chi(g)
    if k is None:
        continue
    value, t = None, time.time()
    for budget in range(k, min(k + 12, g.n) + 1):
        v = run(g, k, budget, rounds=100000, report=10 ** 9)
        if v is True:
            value = budget
            break
    note = ("= k, uniquely colourable" if value == k else
            "= n, vertex-critical" if value == g.n else
            f"= k + {value - k}" if value else "not found in k..k+12")
    print(f"{name:22} {g.n:6} {k:4} {str(value):>5}  {note}"
          f"  [{time.time()-t:.0f}s]", flush=True)
