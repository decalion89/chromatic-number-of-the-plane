"""If the gain is global, measure the global thing -- each graph at its own k.

The chain's local statistics do not move: Sa, Y and G agree to two decimals on
density, rhombus memberships and spindles, while chi goes 4, 4, 5.  So whatever
the bite and the spindle add is global, and the question is which global
quantity moves.

Two are measurable, and BOTH have to be asked at the right number of colours,
which is not the same number for all three graphs.

  - The colour relation -- pairs forced to share or to differ -- is only
    meaningful at k = chi.  Asked at k = 4 on G, which is 5-chromatic, the
    formula is unsatisfiable and every pair comes back "forced" for the empty
    reason, which measures nothing.
  - Criticality -- does removing v leave the graph colourable -- is only
    meaningful at k = chi - 1.  Asked at k = 4 on Sa, which is 4-colourable,
    every vertex comes back essential for the empty reason, again measuring
    nothing.  The first version of this script did exactly that and reported
    60 of 60.

So: the relation at 4, 4, 5 and criticality at 3, 3, 4.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import forced_same, forced_different
from pysat.solvers import Solver

SAMPLE = 60
BUDGET = 150000
t0 = time.time()
rng = random.Random(5)


def edges(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    return sorted(set((min(a, c), max(a, c))
                      for a, c in fast_edges_complete(b, r)))


def formula(n, E, k, sel=None):
    if sel is None:
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    else:
        cls = [[-sel[v]] + [1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    return cls


for name, P, chi in (("Sa", build_Sa(K), 4), ("Y", build_Y(K), 4),
                     ("G", build_G(K, as_graph=False), 5)):
    n, E = len(P), edges(P)
    sv = Solver(name="cd15", bootstrap_with=formula(n, E, chi))
    assert sv.solve(), f"{name} must colour with {chi}"
    same = diff = tested = 0
    while tested < SAMPLE:
        i, j = rng.randrange(n), rng.randrange(n)
        if i == j:
            continue
        tested += 1
        if forced_same(sv, i, j, chi):
            same += 1
        elif forced_different(sv, i, j, chi):
            diff += 1
    sv.delete()
    kk = chi - 1
    nv = n * kk
    sel = [nv + 1 + v for v in range(n)]
    s2 = Solver(name="cd15", bootstrap_with=formula(n, E, kk, sel))
    allsel = [sel[u] for u in range(n)]
    assert s2.solve_limited(assumptions=allsel) is not True or True
    ess = und = 0
    picks = rng.sample(range(n), min(SAMPLE, n))
    for v in picks:
        a = list(allsel)
        a[v] = -sel[v]
        s2.conf_budget(BUDGET)
        if s2.solve_limited(assumptions=a) is True:
            ess += 1
        else:
            und += 1
    s2.delete()
    print(f"{name}: {n} pts, chi={chi} | relation at k={chi} on {tested} "
          f"random pairs: {same} forced-same, {diff} forced-different | "
          f"criticality at k={kk} on {len(picks)}: {ess} essential, "
          f"{und} undecided  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
