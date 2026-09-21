"""Erdos's lattice unioned with a Pythagorean rotation of itself.

Only Y has forced-same pairs here, and Y is a UNION: Sa u rho(Sa).  The base
never forces; the rotation makes it.  So do the same to the lattice, which has
the degree de Grey's family cannot reach.

The rotation is free of any field question if it is Pythagorean.  With
cos = 3/5 and sin = 4/5, rotating the integer grid lands it in (1/5)Z^2, so
scaling everything by five puts the whole union inside Z^2:

    base      5Z^2
    rotated   rho(5i, 5j) = (3i - 4j, 4i + 3j), integer by construction
    edges     squared distance 25r

Every coordinate is an integer, every distance an integer, and the arithmetic
is exact with no field at all.  Within the base the degree is r_2(r) as
before; the cross edges are what the rotation buys.

What is being looked for is a pair forced to AGREE at five colours -- not the
forced-different pairs the ratio has been counting all along.  One of those,
at any distance, spindles to chi >= 6.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from pysat.solvers import Solver
from hn.homcol import agreeing_pairs

SAMPLES = 24
t0 = time.time()


def build(r, half, a=3, b=4, c=5):
    """5Z^2 and its (a,b,c) rotation, inside Z^2, edges at distance^2 = c^2 r."""
    base = [(c * i, c * j) for i in range(-half, half + 1)
            for j in range(-half, half + 1)]
    rot = [(a * i - b * j, b * i + a * j) for i in range(-half, half + 1)
           for j in range(-half, half + 1)]
    pts, idx = [], {}
    for p in base + rot:
        if p not in idx:
            idx[p] = len(pts)
            pts.append(p)
    target = c * c * r
    lim = int(target ** 0.5) + 1
    steps = [(dx, dy) for dx in range(-lim, lim + 1)
             for dy in range(-lim, lim + 1)
             if dx * dx + dy * dy == target and (dx, dy) > (0, 0)]
    E = set()
    for p in pts:
        i = idx[p]
        for dx, dy in steps:
            j = idx.get((p[0] + dx, p[1] + dy))
            if j is not None:
                E.add((min(i, j), max(i, j)))
    return pts, sorted(E), len(base), target


def hunt(r, half, k=5):
    pts, E, nbase, target = build(r, half)
    n = len(pts)
    shared = 2 * nbase - n if n <= 2 * nbase else 0
    Es = set(E)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** r={r}, half={half}: {n} points NOT {k}-COLOURABLE ***",
              flush=True)
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
    print(f"  r={r}, half={half}: {n} points ({shared} shared), {len(E)} "
          f"edges, degree {2*len(E)/n:.2f}; {len(same)} forced-same "
          f"candidates  [{time.time()-t0:.0f}s]", flush=True)
    hits = [(i, j) for i, j in same
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    sv.delete()
    if hits:
        print(f"    *** {len(hits)} FORCED-SAME AT FIVE COLOURS ***",
              flush=True)
        for i, j in hits[:5]:
            dx, dy = pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]
            D = Fr(dx * dx + dy * dy, target)
            print(f"      {pts[i]} and {pts[j]}, squared distance {D} "
                  f"(in units of the edge); spindlable: {4*D >= 1}",
                  flush=True)
    else:
        print(f"    none forced of {len(same)} candidates  "
              f"[{time.time()-t0:.0f}s]", flush=True)


for r, half in ((25, 12), (65, 14), (325, 22)):
    hunt(r, half)
