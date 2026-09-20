"""Are pressure and rho independent?  It matters enormously which.

Sa alone: pressure 3 at four colours, rho 7.  Sa u rot(Sa): rho 5.  If the
pressure there is still 3, then rho fell while the pressure did not move --
and the two are measuring different things.

That would reframe a large part of this branch. Every pressure measurement at
five colours has come back 2, through maximal local enrichment, level-three
densification and unions up to eleven thousand vertices. If pressure and rho
are decoupled, none of that says anything about rho, and a small rho on G is
not excluded by any of it.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure, shrink_forcing_set
from pysat.solvers import Solver

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))


def rho(g, k, rounds=16):
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


pts = list(degrey.build_Sa())
for label in ("Sa", "Sa u rot(Sa)"):
    if label != "Sa":
        seen = set(pts)
        for p in list(pts):
            q = ROT(p)
            if q not in seen:
                seen.add(q)
                pts.append(q)
    g = build_graph(pts)
    deg = g.degrees()
    rel = ColourRelations(g, 4)
    hubs = sorted(range(g.n), key=lambda v: -deg[v])[:4]
    prs = [pressure(rel, v) for v in hubs]
    t = time.time()
    r = rho(g, 4)
    print(f"{label:14}: n={g.n:5}  pressure at the four best pivots {prs}, "
          f"rho = {r}  [{time.time()-t:.0f}s]", flush=True)
print("\nif the pressures match while rho fell, the two are decoupled and every")
print("flat pressure measured at five colours says nothing about rho there.")
