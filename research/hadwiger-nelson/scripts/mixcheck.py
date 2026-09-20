"""The instrument check the ladder depends on: are the samples actually mixed?

A ratio inflates whenever the sampler returns near-identical colourings, and
that is not hypothetical -- one Sa subgraph in forty returns twenty-four
colourings differing on two per cent of vertices and scores 700.  If the
headline figures suffer the same way, the ladder is measuring the solver
rather than the plane.

Two independent proper k-colourings disagree on about (1 - 1/k) of the
vertices.  So measure that, for every graph the ladder rests on, comparing
each pair of samples under the best of the k! colour renamings.
"""
import sys, time, random
from itertools import combinations, permutations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def lattice(side):
    h = F.sqrt(3) / F.rational(2)
    return [Point(F.rational(i) + F.rational(j) / F.rational(2),
                  h * F.rational(j))
            for i in range(side) for j in range(side)]


def check(name, pts, k, pairs_cap=60):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: not colourable", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(7919), []
    for s in range(SAMPLES):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    sv.delete()
    C = np.array(cols, dtype=np.int8)
    nd = 0
    for i in range(n - 1):
        nd += int((C[:, i + 1:] != C[:, i:i + 1]).all(axis=0).sum())
    nd -= len(E)
    ch = n * (n - 1) // 2 * ((k - 1) / k) ** SAMPLES
    perms = list(permutations(range(k)))
    tot, cnt = 0.0, 0
    for a, b in list(combinations(range(SAMPLES), 2))[:pairs_cap]:
        lut = np.empty((len(perms), k), dtype=np.int8)
        for i, p in enumerate(perms):
            lut[i] = p
        best = min(int((C[a] != lut[i][C[b]]).sum()) for i in range(len(perms)))
        tot += best / n
        cnt += 1
    print(f"  {name} at {k}: {n} points, ratio {nd/max(ch,1e-9):.1f}, "
          f"samples disagree on {tot/cnt:.3f} of vertices "
          f"(independent {1-1/k:.3f})  [{time.time()-t0:.0f}s]", flush=True)


check("triangular lattice", lattice(20), 3)
check("Sa", build_Sa(F), 4)
check("Sa", build_Sa(F), 5)
check("Y", build_Y(F), 4)
check("G", build_G(F, as_graph=False), 5)
