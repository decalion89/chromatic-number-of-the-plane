"""Is rho -> k + 1 under unioning a pattern, or a fact about Sa at four?

Sa unioned with a rotated copy drops rho from 7 to 5 at four colours, one
above the floor rho >= k.  If the same happens at three colours on a different
base graph, the pattern is about unioning rather than about Sa, and predicting
rho -> 6 at five colours becomes an extrapolation from two points rather than
a guess from one.

Each base graph is unioned with rotated copies of itself, and rho is measured
the same way throughout: UNSAT cores to get a forcing set fast, then greedy
deletion to make it minimal.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, shrink_forcing_set
from pysat.solvers import Solver


def colourable(g, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


def rho_upper(g, k, rounds=14):
    # rho only means anything when proper k-colourings exist.  Without this
    # guard a union that stops being k-colourable reports every set as
    # "forcing" for want of a counterexample: the triangular patch returned
    # rho = 1 at three colours, below the floor rho >= k, because two copies
    # at the Moser angle are already 4-chromatic -- which is the spindle.
    if not colourable(g, k):
        return "vacuous: not k-colourable"
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


FLD = Field((3, 11))
HALF = FLD.rational(Fraction(1, 2))


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


CASES = [
    ("triangular patch", tri_patch(), 3,
     Rotation(FLD.rational(Fraction(5, 6)),
              FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))),
    ("Moser spindle", moser(), 4,
     Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
              FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))),
]

for name, base, k, rot in CASES:
    print(f"\n{name}, k = {k}  (floor is rho >= {k}):", flush=True)
    cur = list(base)
    for copies in range(1, 6):
        if copies > 1:
            seen, grown = set(cur), list(cur)
            for p in cur:
                q = rot(p)
                if q not in seen:
                    seen.add(q)
                    grown.append(q)
            cur = grown
        g = build_graph(cur)
        t = time.time()
        r = rho_upper(g, k)
        print(f"  {copies} copies: n={g.n:5} m={g.m:6}  rho = "
              f"{r if r is not None else 'not forcing'}"
              + (f"   <- k + 1" if r == k + 1 else "")
              + f"  [{time.time()-t:.0f}s]", flush=True)
        if g.n > 900:
            break
