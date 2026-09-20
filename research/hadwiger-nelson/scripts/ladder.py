"""The decoupling, as a function of the colour count alone, on fixed graphs.

Sa's correlation collapses fifteenfold between four colours and five, with the
graph held fixed.  The obvious next question is whether that is a step or a
slope -- and the triangular lattice answers it at the bottom end, because its
proper 3-colouring is unique up to permuting the colours, so at three colours
every pair of its points is constrained and the ratio should be enormous.

Three graphs, each measured at its own chromatic number and above it, all at
twenty-four samples so the numbers compare:

    lattice patch   chi = 3   at 3, 4, 5
    Sa              chi = 4   at 4, 5, 6
    G               chi = 5   at 5, 6

If the ratio at k = chi falls as chi rises -- unbounded at three, twenty at
four, one at five -- then the difficulty of chi >= 6 is not a matter of
finding the right graph.  It is that the plane stops constraining its own
colourings somewhere between four colours and five.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def lattice(side):
    half = F.rational(1) * F.rational(1) / F.rational(2)
    h = F.sqrt(3) / F.rational(2)
    out = []
    for i in range(side):
        for j in range(side):
            out.append(Point(F.rational(i) + F.rational(j) / F.rational(2),
                             h * F.rational(j)))
    return out


def ratio(name, pts, k):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: NOT {k}-colourable  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(1123), []
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
    pr = n * (n - 1) // 2
    ch = pr * ((k - 1) / k) ** SAMPLES
    print(f"  {name} ({n} pts, {len(E)} edges) at {k} colours: {nd} "
          f"candidates, chance {ch:.1f}, ratio {nd/max(ch,1e-9):.1f}  "
          f"[{time.time()-t0:.0f}s]", flush=True)


L = lattice(20)
for k in (3, 4, 5):
    ratio("triangular lattice", L, k)
Sa = build_Sa(F)
for k in (4, 5, 6):
    ratio("Sa", Sa, k)
G = build_G(F, as_graph=False)
for k in (5, 6):
    ratio("G", G, k)
