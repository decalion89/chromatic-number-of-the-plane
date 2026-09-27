"""Can gadget density be RAISED by construction, or is it pinned?

Everything now turns on one number: critical subgraphs per point.  Sa reaches
0.574 Moser spindles per point and its census collapses to ten; the onset is
about 0.155, measured at two levels.  So the question that decides whether the
programme has any constructive route left is whether that number can be pushed
up deliberately.

Growing for rhombi already pushed rhombus saturation from Sa's 4.47 to 13.65,
three times over -- but rhombi are the raw material, not the gadget.  The
gadget is the SPINDLE, and spindles were never counted along that growth.  So
count them.

Two outcomes, both worth having.  If spindle density climbs well past 0.574,
then gadget density is constructible and the five-colour problem is a question
of doing the same thing with a 509-vertex gadget.  If it saturates near Sa's
own value, that is a structural fact about unit-distance graphs in the plane
and the pinning is not an accident of de Grey's design.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete

t0 = time.time()
rng = random.Random(3)          # the same seed the rhombus growth used
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))


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


def measure(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    assert b.overflow_headroom(r) < 1.0
    n = len(P)
    dm, D2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    nb = [set() for _ in range(n)]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
    rh = sp = 0
    for i in range(n):
        dv = r - r[i]
        s = b._field_square(dv[:, :dm]) + b._field_square(dv[:, dm:])
        g = s[:, 0] == 3 * D2
        for m in range(1, dm):
            g &= s[:, m] == 0
        far = [int(j) for j in np.nonzero(g)[0]]
        for j in far:
            if j > i:
                sh = len(nb[i] & nb[j])
                if sh >= 2:
                    rh += sh * (sh - 1) // 2
        for a in range(len(far) - 1):
            for c in range(a + 1, len(far)):
                x, y = far[a], far[c]
                if y not in nb[x]:
                    continue
                sx, sy = sorted(nb[i] & nb[x]), sorted(nb[i] & nb[y])
                if len(sx) < 2 or len(sy) < 2:
                    continue
                if any(len({i, x, y, sx[p], sx[q], sy[u], sy[v]}) == 7
                       for p in range(len(sx) - 1)
                       for q in range(p + 1, len(sx))
                       for u in range(len(sy) - 1)
                       for v in range(u + 1, len(sy))):
                    sp += 1
    return n, len(E), rh, sp, nb, b, r, E


P = build_Sa(K)
n, ne, rh, sp, nb, b, r, E = measure(P)
print(f"start: {n} pts, {ne} edges, {rh} rhombi, {sp} spindles "
      f"({sp/n:.3f}/pt)  [{time.time()-t0:.0f}s]", flush=True)
BATCH = 40
for step in range(1, 16):
    cand = candidates(P, 30000)
    if not cand:
        break
    gb = IntBasis.covering(P + cand)
    dim, D2 = gb.dim, gb.D * gb.D
    Pr = gb.rows(P)
    scored = []
    for q in cand:
        o = gb.rows([q])[0]
        d = Pr - o
        sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
        hit = sq[:, 0] == D2
        for m in range(1, dim):
            hit &= sq[:, m] == 0
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
        far = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
        g3 = far[:, 0] == 3 * D2
        for m in range(1, dim):
            g3 &= far[:, m] == 0
        for w in np.nonzero(g3)[0]:
            sh = len(nbq & nb[int(w)])
            gain += sh * (sh - 1) // 2
        scored.append((gain, len(nbq), q))
    scored.sort(key=lambda t: (-t[0], -t[1]))
    add = [q for g, dgr, q in scored[:BATCH] if g > 0]
    if not add:
        break
    P = list(dict.fromkeys(P + add))
    n, ne, rh, sp, nb, b, r, E = measure(P)
    print(f"  step {step}: {n} pts, {ne} edges ({ne/n:.2f}/v), {rh} rhombi "
          f"({rh*4/n:.2f} memb/pt), {sp} spindles ({sp/n:.3f}/pt)"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
