"""Does rho actually DROP when the graph grows, or only fail to rise?

Monotonicity gives rho(W) <= rho(G) for G inside W, which is why the strategy
turned towards bigger graphs.  But <= is not <, and the whole plan rests on
the difference: if enlarging a graph leaves rho where it was, there is no
route from 1581 down to 63 however much is added.

Measured where rho is cheap: Sa at FOUR colours, where the core method runs in
under a second.  Sa is enlarged by rotating it about the origin through the
half-Moser angle and taking unions, so each graph strictly contains the last
and the chromatic number stays 4.

The core method returns a forcing set, not necessarily the smallest, so the
numbers are upper bounds on rho -- but they are computed the same way at every
size, so the TREND is the thing, and a trend that flattens is the answer.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))


def shrink(graph, k, rounds=14):
    n = graph.n
    NV = n * k

    def x(v, c):
        return 1 + v * k + c

    def sel(v):
        return NV + 1 + v

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in graph.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    for v in range(n):
        cls.append([-sel(v), -x(v, 0)])

    S = list(range(n))
    for _ in range(rounds):
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve(assumptions=[sel(v) for v in S]):
                return S, False
            core = sorted({abs(l) - NV - 1 for l in s.get_core()})
        if len(core) >= len(S):
            return core, True
        S = core
    return S, True


pts = list(degrey.build_Sa())
print(f"{'copies':>7} {'n':>6} {'m':>7} {'minimal rho':>12}   time")
print("-" * 44)
cur = list(pts)
for copies in range(1, 8):
    if copies > 1:
        seen = set(cur)
        grown = list(cur)
        for p in cur:
            q = ROT(p)
            if q not in seen:
                seen.add(q)
                grown.append(q)
        cur = grown
    g = build_graph(cur)
    t = time.time()
    S, ok = shrink(g, 4)
    # the core is forcing but not minimal, and the plateau could be an
    # artefact of that -- so peel it down to a MINIMAL forcing set, which is
    # cheap once the set is small
    if ok:
        from hn.forced import ColourRelations, shrink_forcing_set
        S = shrink_forcing_set(ColourRelations(g, 4), S)
    print(f"{copies:>7} {g.n:>6} {g.m:>7} {len(S):>12}   "
          f"[{time.time()-t:.0f}s]" + ("" if ok else "  (not forcing)"),
          flush=True)
