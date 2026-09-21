"""Grow for saturation, which is the one variable measured to move the census.

Everything else grown here was grown for edges, and edges do not move the
quantity of interest: density rises by a factor of two across the gap between
criticality and rigidity while rhombus saturation rises by eight, and choosing
points for saturation moves the census off a ceiling that choosing them at
random never leaves.

So grow for it directly.  Candidates are the points at unit distance from two
points already present -- where unit circles meet, and the only place a new
point can attach to more than one vertex.  Each is scored by the RHOMBI it
completes, not the edges it adds, and the best is taken.

Two things are watched.  Saturation, to see how far past Sa's 4.47 this can
go; and 4-colourability, because a saturation-grown graph that stops colouring
below 1581 points is a 5-chromatic unit-distance graph smaller than G.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
rng = random.Random(3)
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


def analyse(P):
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
    per = defaultdict(int)
    rh = 0
    for i in range(n):
        dv = r - r[i]
        s = b._field_square(dv[:, :dm]) + b._field_square(dv[:, dm:])
        good = s[:, 0] == 3 * D2
        for m in range(1, dm):
            good &= s[:, m] == 0
        for j in np.nonzero(good)[0]:
            j = int(j)
            if j <= i:
                continue
            sh = nb[i] & nb[j]
            if len(sh) >= 2:
                c2 = len(sh) * (len(sh) - 1) // 2
                rh += c2
                per[i] += c2
                per[j] += c2
                for x in sh:
                    per[x] += len(sh) - 1
    return n, len(E), rh, sum(per.values()) / n, b, r, E


def colours(P, E, k):
    n = len(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    out = s.solve()
    s.delete()
    return out


P = build_Sa(K)
n, ne, rh, sat, b, r, E = analyse(P)
print(f"start from Sa: {n} pts, {ne} edges, {rh} rhombi, {sat:.2f} per point"
      f"  [{time.time()-t0:.0f}s]", flush=True)
BATCH = 40
for step in range(1, 26):
    cand = candidates(P, 30000)
    if not cand:
        break
    gb = IntBasis.covering(P + cand)
    dim, D2 = gb.dim, gb.D * gb.D
    Pr = gb.rows(P)
    nb = [set() for _ in range(len(P))]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
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
        # A rhombus is a pair at squared distance 3 with two common unit
        # neighbours.  A new point q completes one in two distinct ways, and
        # both have to be counted:
        #   (a) q is one of the two COMMON NEIGHBOURS -- some pair inside q's
        #       own neighbourhood is at squared distance 3, and each common
        #       neighbour it already has pairs with q to make a new rhombus;
        #   (b) q is one END of the pair -- some existing point is at squared
        #       distance 3 from q, and any two of their common neighbours
        #       close a rhombus.
        gain = 0
        nl = sorted(nbq)
        for ii in range(len(nl) - 1):
            di = Pr[nl[ii + 1:]] - Pr[nl[ii]]
            s3 = gb._field_square(di[:, :dim]) + gb._field_square(di[:, dim:])
            g3 = s3[:, 0] == 3 * D2
            for m in range(1, dim):
                g3 &= s3[:, m] == 0
            for off in np.nonzero(g3)[0]:
                jj = nl[ii + 1 + int(off)]
                gain += len(nb[nl[ii]] & nb[jj])            # (a)
        far = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
        g3 = far[:, 0] == 3 * D2
        for m in range(1, dim):
            g3 &= far[:, m] == 0
        for w in np.nonzero(g3)[0]:
            sh = len(nbq & nb[int(w)])
            gain += sh * (sh - 1) // 2                      # (b)
        scored.append((gain, len(nbq), q))
    scored.sort(key=lambda t: (-t[0], -t[1]))
    add = [q for g, dgr, q in scored[:BATCH] if g > 0]
    if not add:
        print("  no candidate completes a rhombus, stop", flush=True)
        break
    P = list(dict.fromkeys(P + add))
    n, ne, rh, sat, b, r, E = analyse(P)
    ok4 = colours(P, E, 4)
    print(f"  step {step}: {n} pts, {ne} edges ({ne/n:.2f}/v), {rh} rhombi, "
          f"{sat:.2f} per point, 4-colourable={ok4}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok4:
        print(f"*** {n} points NOT 4-COLOURABLE -- 5-chromatic, against "
              f"G's 1581 ***", flush=True)
        break
print("DONE", flush=True)
