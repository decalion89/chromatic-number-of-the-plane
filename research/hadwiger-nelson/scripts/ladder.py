"""Is the ladder real?  A uniquely 3-colourable graph gave a 4-chromatic one.

The triangular lattice is uniquely 3-colourable -- rho = 3 -- and unioning it
with a Moser-rotated copy makes it 4-chromatic, measured at three radii.  That
is one rung.

The rung above needs the union to be uniquely 4-COLOURABLE, rho = 4.  If it
is, a third copy at the right angle gives chi >= 5 by the same move, and the
step after that gives 6.  If it is not, the ladder stops at the first rung and
the reason is worth having.

So: rho of the union, by decision, budget swept up from the floor.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fractions import Fraction
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
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


for r in (2, 3):
    base = tri_patch(r)
    seen, both = set(base), list(base)
    for p in base:
        q = MOSER(p)
        if q not in seen:
            seen.add(q)
            both.append(q)
    w = build_graph(both)
    k = 4 if colourable(w, 4) and not colourable(w, 3) else None
    print(f"\npatch radius {r} + Moser copy: n={w.n}, m={w.m}, chi = {k}",
          flush=True)
    if k != 4:
        continue
    for budget in (4, 5, 6, 7, 8):
        t = time.time()
        v = run(w, 4, budget, rounds=200000, report=10 ** 9)
        print(f"  budget {budget}: forcing set exists = {v}"
              + ("   *** rho = 4, UNIQUELY 4-COLOURABLE ***"
                 if v is True and budget == 4 else "")
              + f"  [{time.time()-t:.0f}s]", flush=True)
        if v is True:
            print(f"  => rho = {budget}"
                  + ("  the ladder continues" if budget == 4
                     else f", which is k + {budget - 4}"), flush=True)
            break
