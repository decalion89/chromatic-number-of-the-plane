"""Can unioning reach the floor rho = k, and does the drop happen on other bases?

Sa u rot(Sa) reaches rho = 5 at four colours, one above the floor.  Two
questions follow.  Does the same happen to de Grey's other 4-chromatic pieces,
Sb and Y -- or is it a fact about Sa?  And can a more aggressive union reach
rho = 4 exactly, which would mean a uniquely 4-colourable unit-distance graph?

Every measurement is guarded: a union that stops being 4-colourable makes rho
vacuous, and reports as such rather than as a triumph.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, shrink_forcing_set
from pysat.solvers import Solver

FLD = degrey.DEGREY_FIELD
HALF = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
                FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))
FULL = Rotation(FLD.rational(Fraction(5, 6)),
                FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
SIXTY = Rotation(FLD.rational(Fraction(1, 2)),
                 FLD.sqrt(3) * FLD.rational(Fraction(1, 2)))


def colourable(g, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


def rho(g, k, rounds=16):
    if not colourable(g, k):
        return None
    n, NV = g.n, g.n * k

    def x(v, c):
        return 1 + v * k + c

    def sel(v):
        return NV + 1 + v

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    for v in range(n):
        cls.append([-sel(v), -x(v, 0)])
    S = list(range(n))
    for _ in range(rounds):
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve(assumptions=[sel(v) for v in S]):
                return None
            core = sorted({abs(l) - NV - 1 for l in s.get_core()})
        if len(core) >= len(S):
            break
        S = core
    return len(shrink_forcing_set(ColourRelations(g, k), S))


def grow(pts, rots):
    seen, out = set(pts), list(pts)
    for r in rots:
        for p in list(pts):
            q = r(p)
            if q not in seen:
                seen.add(q)
                out.append(q)
    return out


print("de Grey's 4-chromatic pieces, each with one rotated copy:", flush=True)
for name, base in (("Sa", degrey.build_Sa()), ("Sb", degrey.build_Sb()),
                   ("Y", degrey.build_Y())):
    for label, rots in (("alone", []), ("+half", [HALF])):
        pts = grow(list(base), rots)
        g = build_graph(pts)
        t = time.time()
        r = rho(g, 4)
        print(f"  {name:3} {label:6}: n={g.n:5} m={g.m:6}  rho = "
              f"{r if r is not None else 'vacuous / not forcing'}"
              f"  [{time.time()-t:.0f}s]", flush=True)

print("\npushing Sa harder, to see whether the floor rho = 4 is reachable:",
      flush=True)
for label, rots in (("+half", [HALF]), ("+half,full", [HALF, FULL]),
                    ("+half,full,60", [HALF, FULL, SIXTY]),
                    ("+half,full,60,half^2", [HALF, FULL, SIXTY,
                                              lambda p: HALF(HALF(p))])):
    pts = grow(list(degrey.build_Sa()), rots)
    g = build_graph(pts)
    t = time.time()
    r = rho(g, 4)
    print(f"  Sa {label:22}: n={g.n:5} m={g.m:6}  rho = "
          f"{r if r is not None else 'vacuous / not forcing'}"
          + ("   *** FLOOR ***" if r == 4 else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
