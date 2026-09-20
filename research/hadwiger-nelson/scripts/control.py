"""The control the conflict jump needs.

A translate union of G costs 723211 conflicts to sweep against G's own 6410,
and its dearest pair 3689 against 29.  Before reading that as tightening, rule
out the dull explanation: the union has 3026 vertices against 1581, and a
bigger instance is harder to solve even when it is just as loose.

So take two copies of G that cannot possibly be tightening each other -- the
same point set translated so far away that no cross edge is geometrically
possible -- and sweep that.  If it also costs hundreds of thousands, the jump
is size.  If it costs about twice G's own, the jump is real.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()
pts = build_G(F, as_graph=False)
GE = list(build_graph(pts).edges())
zf = [(float(p.x), float(p.y)) for p in pts]
n0 = len(pts)
SHIFT = 1000.0
allz = zf + [(x + SHIFT, y) for x, y in zf]
E = [(a, b) for a, b in GE] + [(a + n0, b + n0) for a, b in GE]
n = 2 * n0
print(f"control: {n} points, {len(E)} edges, no cross edge possible "
      f"(shifted by {SHIFT})  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * 5 + c for c in range(5)] for v in range(n)]
for a, b in E:
    for c in range(5):
        cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
sv.solve()
cand = []
for i in range(n):
    ai, bi = allz[i]
    for j in range(i + 1, n):
        v = (ai - allz[j][0]) ** 2 + (bi - allz[j][1]) ** 2
        if v > 36.0:
            continue
        D = Fr(round(v * 1584), 1584)
        if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
            continue
        cand.append((i, j))


def conf():
    return sv.accum_stats().get("conflicts", 0)


base, dear, hard = conf(), 0, 0
for i, j in cand:
    sv.conf_budget(60000)
    r = sv.solve_limited(assumptions=[1 + i * 5, -(1 + j * 5)])
    if r is None:
        hard += 1
    c2 = conf()
    dear = max(dear, c2 - base)
    base = c2
print(f"  {len(cand)} pairs, {conf()} conflicts, dearest {dear}, {hard} "
      f"over budget  [{time.time()-t0:.0f}s]", flush=True)
print(f"  compare: G alone 6410 and 29; the translate union 723211 and 3689",
      flush=True)
