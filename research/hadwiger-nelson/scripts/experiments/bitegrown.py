"""Bite the grown graph, where the bite has more ring to grab.

Raising local density did nothing to the census -- ten of ten at every step of
a growth that quadrupled the spindle count.  But the bite is not a local
operation, and what it needs is not density: it needs a RING, and it takes the
census from ten to three on Sa's six-point one.

The grown graph has the same origin and the same rings, thickened.  So the
question the growth was worth asking after all is not whether density tightens
anything by itself -- it does not -- but whether it gives the bite more to
work with.  If Sa u rho(Sa) reads three, does grown u rho(grown) read less?

Both are 4-colourable, so a reading of zero is the graph ceasing to colour: a
5-chromatic unit-distance graph at about two thousand points by a route that
is not de Grey's spindle.  Measured against Sa's three at every ring the field
can join about the origin.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 4
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
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
dm, d2 = b0.dim, b0.D * b0.D
zi = P0.index(Point(K.zero(), K.zero()))
d = r0 - r0[zi]
sq = b0._field_square(d[:, :dm]) + b0._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring0 = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
ring0.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
SEVEN = [P0[zi]] + [P0[i] for i in ring0]
ZERO = P0[zi]

P = list(P0)
for step in range(STEPS):
    cand = candidates(P, 30000)
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
        scored.append((gain, len(nbq), q))
    scored.sort(key=lambda t: (-t[0], -t[1]))
    P = list(dict.fromkeys(P + [q for g, dg, q in scored[:40] if g > 0]))
print(f"grown to {len(P)} points in {STEPS} steps  [{time.time()-t0:.0f}s]",
      flush=True)

Ds = [D for D in sorted({Fr(x, y) for y in range(1, 10)
                         for x in range(1, 10 * y + 1)} - {Fr(1)})
      if closable_distance(D)]
print(f"{len(Ds)} rings about the origin  [{time.time()-t0:.0f}s]", flush=True)


def census(U, label):
    b = IntBasis.covering(U)
    r = b.rows(U)
    if b.overflow_headroom(r) >= 1.0:
        return None
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(U)
    W = [U.index(q) for q in SEVEN]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        print(f"*** {label}: {n} points NOT 4-COLOURABLE ***", flush=True)
        return 0
    keep = sum(1 for p in TEN
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(p)
                                        for e in blk]))
    sv.delete()
    return keep


print(f"unbitten grown graph: census {census(list(P), 'grown')} of 10",
      flush=True)
best = 10
for D in Ds:
    rot = rotation_joining(D, K)
    for sgn in (+1, -1):
        f = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(ZERO)
        s2 = set(P)
        U = list(P) + [q for q in (f(p) for p in P) if q not in s2]
        c = census(U, f"grown bitten at D={D} dir {sgn}")
        if c is None:
            continue
        if c < best:
            best = c
            print(f"*** D={D}, dir {sgn}: {len(U)} points, census {c} of 10 "
                  f"-- new best (Sa's bite gives 3) ***"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
        if c == 0:
            sys.exit()
print(f"\nbest census {best} of 10  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
