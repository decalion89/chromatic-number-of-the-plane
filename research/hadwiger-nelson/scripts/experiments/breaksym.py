"""Break the symmetry and see how far below three the census goes.

The three surviving patterns of Sa u rho(Sa) form one orbit under sixty-degree
rotation, so any union that keeps that symmetry has three or none.  Adding an
ASYMMETRIC point cannot add patterns -- adding vertices never does -- but it
can kill part of the orbit, and a census of one would mean the centre and the
three antipodal pairs are pinned up to renaming colours: strictly more forced
structure than de Grey's construction uses.

Candidates are the points at unit distance from two points already present,
which is where unit circles meet and the only place a new point can attach to
more than one vertex.  A candidate is screened by testing the three patterns:
all three satisfiable is the fast answer and means nothing happened.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Sb
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
t0 = time.time()
rng = random.Random(7)
THREE = [
    [[0, 1, 4], [2, 3, 5, 6]],
    [[1, 2, 4, 5], [0, 3, 6]],
    [[0, 2, 5], [1, 3, 4, 6]],
]
ZERO = Point(K.zero(), K.zero())
U, seen = [], set()
for p in build_Sa(K) + build_Sb(K):
    if p not in seen:
        seen.add(p)
        U.append(p)
print(f"Sa u Sb: {len(U)} points  [{time.time()-t0:.0f}s]", flush=True)

CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]


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


half = K.rational(Fr(1, 2))
cands, cs = [], set(U)
for _ in range(40000):
    A, B = U[rng.randrange(len(U))], U[rng.randrange(len(U))]
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
        if q not in cs:
            cs.add(q)
            cands.append(q)
print(f"{len(cands)} candidate points  [{time.time()-t0:.0f}s]", flush=True)


def seven(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    ci = P.index(ZERO)
    key = {tuple(r[i]): i for i in range(len(P))}
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
    ms = set(ring)
    pairs = [(j, key[tuple(2 * r[ci] - r[j])]) for j in ring
             if tuple(2 * r[ci] - r[j]) in key
             and key[tuple(2 * r[ci] - r[j])] in ms
             and key[tuple(2 * r[ci] - r[j])] > j]
    six = [x for pr in pairs[:3] for x in pr]
    six.sort(key=lambda i: math.atan2(float(P[i].y), float(P[i].x)))
    return len(P), E, [ci] + six


# how many neighbours does each candidate have?  more contact first
b0 = IntBasis.covering(U)
r0 = b0.rows(U)
scored = []
gb = IntBasis.covering(U + cands[:3000])
dim, D2 = gb.dim, gb.D * gb.D
Ur = gb.rows(U)
for q in cands[:3000]:
    o = gb.rows([q])[0]
    d = Ur - o
    sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
    hit = sq[:, 0] == D2
    for m in range(1, dim):
        hit &= sq[:, m] == 0
    scored.append((int(hit.sum()), q))
scored.sort(key=lambda t: -t[0])
print(f"best candidate has {scored[0][0]} neighbours; trying the top ones"
      f"  [{time.time()-t0:.0f}s]", flush=True)

for deg, q in scored[:25]:
    V = U + [q]
    n, E, W = seven(V)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"*** adding a {deg}-neighbour point makes {n} points NOT "
              f"4-COLOURABLE ***  [{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        break
    keep = sum(1 for part in THREE
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(part)
                                        for e in blk]))
    sv.delete()
    tag = " *** BELOW THREE ***" if keep < 3 else ""
    print(f"  +1 point (degree {deg}): {n} pts, census {keep} of 3{tag}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
