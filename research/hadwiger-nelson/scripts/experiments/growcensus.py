"""Does raising gadget density past Sa's tighten the census further?

Growing for rhombi raises the spindle count per point from Sa's 0.574 to 1.86
by step eight -- three times over, and well past the 0.155 onset and the 0.574
collapse.  If gadget density is what rigidity tracks, the census should keep
falling.  The grown graph CONTAINS Sa, so its census on Sa's own seven points
is at most ten already; the question is whether it goes below, and how far.

It cannot reach zero without the graph ceasing to 4-colour, and the growth was
measured to 4-colour at every step, so one pattern at least survives.  Between
ten and one is the whole of the remaining signal, and it is ten SAT calls to
read.

Regrown with the same seed, so the trajectory matches the density measurement
step for step.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
rng = random.Random(3)
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))
TEN = [
    [[0, 1, 2, 3, 4, 5, 6]],
    [[1, 2], [0, 3, 4, 5, 6]], [[2, 3], [0, 1, 4, 5, 6]],
    [[0, 1, 2, 3, 4], [5, 6]], [[3, 4], [0, 1, 2, 5, 6]],
    [[0, 1, 4], [2, 3, 5, 6]], [[0, 2, 3, 4, 5], [1, 6]],
    [[1, 2, 4, 5], [0, 3, 6]], [[4, 5], [0, 1, 2, 3, 6]],
    [[0, 2, 5], [1, 3, 4, 6]],
]


def rsqrt(e):
    if any(x for x in e.c[1:]):
        return None
    v = e.c[0]
    if v <= 0:
        return None
    for rr in CLASSES:
        q = v / rr
        nn, dd = q.numerator, q.denominator
        rn, rd = int(round(nn ** .5)), int(round(dd ** .5))
        if rn * rn == nn and rd * rd == dd:
            s = K.rational(Fr(rn, rd))
            return s if rr == 1 else K.sqrt(rr) * s
    return None


def candidates(P, tries):
    out, have = [], set(P)
    n = len(P)
    for _ in range(tries):
        A, B = P[rng.randrange(n)], P[rng.randrange(n)]
        if A == B:
            continue
        D = A.dist2(B)
        if not .05 < float(D) < 3.99:
            continue
        sD, s4 = rsqrt(D), rsqrt(K.rational(4) - D)
        if sD is None or s4 is None:
            continue
        inv = K.rational(1) / sD
        mx, my = (A.x + B.x) * half, (A.y + B.y) * half
        nx = -(B.y - A.y) * inv * s4 * half
        ny = (B.x - A.x) * inv * s4 * half
        for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if q not in have:
                have.add(q)
                out.append(q)
    return out


P0 = build_Sa(K)
b0 = IntBasis.covering(P0)
r0 = b0.rows(P0)
dm0, d20 = b0.dim, b0.D * b0.D
zi = P0.index(Point(K.zero(), K.zero()))
d = r0 - r0[zi]
sq = b0._field_square(d[:, :dm0]) + b0._field_square(d[:, dm0:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm0):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d20) == 4]
ring.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
SEVEN = [P0[zi]] + [P0[i] for i in ring]
k = 4


def census(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    W = [P.index(q) for q in SEVEN]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None, len(E)
    keep = [p for p in TEN
            if sv.solve(assumptions=[1 + W[e] * k + bi
                                     for bi, blk in enumerate(p)
                                     for e in blk])]
    sv.delete()
    return keep, len(E)


P = list(P0)
keep, ne = census(P)
print(f"Sa: {len(P)} pts, {ne} edges, census {len(keep)} of 10"
      f"  [{time.time()-t0:.0f}s]", flush=True)
BATCH = 40
nb_cache = None
for step in range(1, 13):
    cand = candidates(P, 30000)
    if not cand:
        break
    gb = IntBasis.covering(P + cand)
    dim, D2 = gb.dim, gb.D * gb.D
    Pr = gb.rows(P)
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(gb, Pr)))
    nb = [set() for _ in range(len(P))]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
    scored = []
    for q in cand:
        o = gb.rows([q])[0]
        dd = Pr - o
        s = gb._field_square(dd[:, :dim]) + gb._field_square(dd[:, dim:])
        hit = s[:, 0] == D2
        for m in range(1, dim):
            hit &= s[:, m] == 0
        nbq = set(int(x) for x in np.nonzero(hit)[0])
        if len(nbq) < 2:
            continue
        gain = 0
        nl = sorted(nbq)
        for ii in range(len(nl) - 1):
            di = Pr[nl[ii + 1:]] - Pr[nl[ii]]
            s3 = gb._field_square(di[:, :dim]) + gb._field_square(di[:, dim:])
            g3 = s3[:, 0] == 3 * D2
            for m in range(1, dim):
                g3 &= s3[:, m] == 0
            for off in np.nonzero(g3)[0]:
                gain += len(nb[nl[ii]] & nb[nl[ii + 1 + int(off)]])
        g3 = s[:, 0] == 3 * D2
        for m in range(1, dim):
            g3 &= s[:, m] == 0
        for w in np.nonzero(g3)[0]:
            sh = len(nbq & nb[int(w)])
            gain += sh * (sh - 1) // 2
        scored.append((gain, len(nbq), q))
    scored.sort(key=lambda t: (-t[0], -t[1]))
    add = [q for g, dgr, q in scored[:BATCH] if g > 0]
    if not add:
        break
    P = list(dict.fromkeys(P + add))
    keep, ne = census(P)
    if keep is None:
        print(f"*** step {step}: {len(P)} pts NOT 4-COLOURABLE ***",
              flush=True)
        break
    print(f"  step {step}: {len(P)} pts, {ne} edges, census {len(keep)} of 10"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if len(keep) <= 3:
        print(f"      survivors: {keep}", flush=True)
print("DONE", flush=True)
