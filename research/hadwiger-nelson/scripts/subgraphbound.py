"""rho <= |smallest k-chromatic subgraph|, and where the union beats it.

Two routes to a small rho, and they are different in kind.

STRUCTURAL. A k-chromatic subgraph uses all k colours in every k-colouring of
the whole graph, so rho <= its size.  At four colours the smallest
unit-distance one is the Moser spindle, seven vertices -- and rho(Sa,4) is
exactly 7, so Sa sits right on that bound even though the minimal forcing set
it actually admits is a different, 3-chromatic, nearly edgeless five... no,
seven.

AMBIENT. Unioning Sa with a rotated copy takes rho to 5, BELOW the structural
bound, with a set that contains no 4-chromatic subgraph at all. That is the
interesting mechanism, and the only one that could ever reach 63 at five
colours, where the structural bound is the smallest 5-chromatic unit-distance
graph -- around five hundred vertices in the published record.

This measures how far below the structural bound the ambient mechanism gets,
at four colours, where both numbers are computable.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, shrink_forcing_set
from pysat.solvers import Solver

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))


def smallest_kchromatic_subgraph(g, k, tries=200):
    """A k-chromatic subgraph, shrunk greedily by an UNSAT core."""
    n, NV = g.n, g.n * (k - 1)

    def x(v, c):
        return 1 + v * (k - 1) + c

    def sel(v):
        return NV + 1 + v

    cls = [[x(v, c) for c in range(k - 1)] for v in range(n)]
    for a, b in g.edges():
        for c in range(k - 1):
            cls.append([-x(a, c), -x(b, c)])
    # selector switches a vertex's "must be coloured" clause on
    cls = [[-sel(v)] + [x(v, c) for c in range(k - 1)] for v in range(n)] + \
          [cl for cl in cls if len(cl) == 2]
    S = list(range(n))
    for _ in range(12):
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve(assumptions=[sel(v) for v in S]):
                return None
            core = sorted({abs(l) - NV - 1 for l in s.get_core()})
        if len(core) >= len(S):
            break
        S = core
    return len(S)


def rho(g, k):
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
    for _ in range(16):
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
    t = time.time()
    sub = smallest_kchromatic_subgraph(g, 4)
    r = rho(g, 4)
    print(f"{label:14}: n={g.n:5}  smallest 4-chromatic subgraph found "
          f"{sub},  rho = {r}"
          + ("   <- ambient beats structural" if r is not None and sub
             and r < sub else "")
          + f"  [{time.time()-t:.0f}s]", flush=True)
