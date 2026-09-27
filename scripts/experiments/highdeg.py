"""Escape the degree ceiling: Erdos's lattice, where degree is arithmetic.

Average degree tops out near eleven in de Grey's family -- 9.96 at G, 10.64 at
G*, 11.36 at Z -- and the four-colour curve says correlation rises steeply
with it while the five-colour curve is flat up to that ceiling.  The ceiling
is the family's, not the plane's.

Erdos's construction lifts it.  Take the integer grid and scale by 1/sqrt(r):
two points are at distance exactly 1 when their difference satisfies
dx^2 + dy^2 = r, so the degree of an interior point is r_2(r), the number of
representations of r as a sum of two squares -- which is 4 times the product
of (a_i + 1) over the primes 1 mod 4 dividing r.  Choose r and the degree
follows:

    r =   25 = 5^2            12
    r =   65 = 5.13           16
    r =  325 = 5^2.13         24
    r = 1105 = 5.13.17        32

Everything is exact in integers, with no field at all, and the graph is
genuinely unit-distance after the scaling.  Nor does it need to be 5-chromatic
for the measurement: the correlation ratio is defined for any graph that
colours.

If degree is the lever, a graph at sixteen or twenty-four should show at five
colours what nothing at eleven has shown.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from itertools import combinations, permutations
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def erdos(r, side):
    """The grid of the given side, edges at squared distance r."""
    pts = [(i, j) for i in range(side) for j in range(side)]
    idx = {p: i for i, p in enumerate(pts)}
    steps = []
    lim = int(r ** 0.5) + 1
    for dx in range(-lim, lim + 1):
        for dy in range(-lim, lim + 1):
            if dx * dx + dy * dy == r and (dx, dy) > (0, 0):
                steps.append((dx, dy))
    E = set()
    for (x, y) in pts:
        for dx, dy in steps:
            q = (x + dx, y + dy)
            j = idx.get(q)
            if j is not None:
                a, b = idx[(x, y)], j
                E.add((min(a, b), max(a, b)))
    return pts, sorted(E)


def measure(name, pts, E, k):
    n = len(pts)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: NOT {k}-COLOURABLE  "
              f"[{time.time()-t0:.0f}s]", flush=True)
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
    lut = np.array(list(permutations(range(k))), dtype=np.int8)
    tot, cnt = 0.0, 0
    for a, b in list(combinations(range(SAMPLES), 2))[:30]:
        tot += min(int((C[a] != lut[i][C[b]]).sum())
                   for i in range(len(lut))) / n
        cnt += 1
    print(f"  {name} at {k}: {n} points, {len(E)} edges, degree "
          f"{2*len(E)/n:.2f}, ratio {nd/max(ch,1e-9):.1f}, spread "
          f"{tot/cnt:.3f} (independent {1-1/k:.3f})  "
          f"[{time.time()-t0:.0f}s]", flush=True)


for r, side in ((25, 30), (65, 40), (325, 80)):
    pts, E = erdos(r, side)
    print(f"r={r}, grid {side}x{side}: {len(pts)} points, {len(E)} edges, "
          f"degree {2*len(E)/len(pts):.2f}  [{time.time()-t0:.0f}s]",
          flush=True)
    for k in (4, 5):
        measure(f"Erdos r={r}", pts, E, k)
