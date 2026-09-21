"""Forced pairs in Erdos's lattice at five colours, where the correlation is.

The lattice scaled by 1/sqrt(r) has degree r_2(r) and, measured, a correlation
ratio at FIVE colours of 8.7 at degree 9.57, 12.7 at 12.18 and 26.4 at 17.64
-- where de Grey's family gives 1.1 at degree 9.96.  Twenty-six is the level
at which forcing appeared one level down.

And there is no field obstruction here.  Every squared distance is m/r with m
an integer, so it is rational; the spindle of a forced pair at squared
distance D exists in the reals for any D >= 1/4, and the rotated copy has
coordinates in Q(sqrt r, sqrt(4D-1)), which is a perfectly good real field.
Nothing has to lie in anyone's chain.

So: sample, filter, and put every survivor to the solver.  A forced pair at
five colours, at any distance at all, spindles to chi >= 6.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from pysat.solvers import Solver
from hn.homcol import agreeing_pairs

SAMPLES = 24
t0 = time.time()


def erdos(r, side):
    pts = [(i, j) for i in range(side) for j in range(side)]
    idx = {p: i for i, p in enumerate(pts)}
    lim = int(r ** 0.5) + 1
    steps = [(dx, dy) for dx in range(-lim, lim + 1)
             for dy in range(-lim, lim + 1)
             if dx * dx + dy * dy == r and (dx, dy) > (0, 0)]
    E = set()
    for (x, y) in pts:
        a = idx[(x, y)]
        for dx, dy in steps:
            j = idx.get((x + dx, y + dy))
            if j is not None:
                E.add((min(a, j), max(a, j)))
    return pts, sorted(E)


def hunt(r, side, k=5):
    pts, E = erdos(r, side)
    n = len(pts)
    Es = set(E)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  r={r}: NOT {k}-COLOURABLE on {n} points -- that is the "
              f"answer  [{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(4242), []
    for s in range(SAMPLES):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    same = [(i, j) for i, j in agreeing_pairs(cols, colours=k)
            if (i, j) not in Es]
    C = np.array(cols, dtype=np.int8)
    diff = []
    for i in range(n - 1):
        msk = (C[:, i + 1:] != C[:, i:i + 1]).all(axis=0)
        for off in np.nonzero(msk)[0]:
            j = int(off) + i + 1
            if (i, j) not in Es:
                diff.append((i, j))
    print(f"  r={r}: {n} points, {len(E)} edges, degree {2*len(E)/n:.2f}; "
          f"{len(same)} forced-same candidates, {len(diff)} forced-different "
          f"candidates  [{time.time()-t0:.0f}s]", flush=True)
    fs = [(i, j) for i, j in same
          if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    fd = [(i, j) for i, j in diff[:3000]
          if not sv.solve(assumptions=[1 + i * k, 1 + j * k])]
    sv.delete()
    print(f"    VERIFIED at five colours: {len(fs)} forced-same, {len(fd)} "
          f"forced-different of {min(len(diff),3000)} checked  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for i, j in fs[:5]:
        dx = pts[i][0] - pts[j][0]
        dy = pts[i][1] - pts[j][1]
        D = Fr(dx * dx + dy * dy, r)
        print(f"    *** FORCED SAME: {pts[i]} and {pts[j]}, squared "
              f"distance {D} -- spindlable iff 4D >= 1: {4*D >= 1} ***",
              flush=True)


for r, side in ((25, 30), (65, 40), (325, 80)):
    hunt(r, side)
