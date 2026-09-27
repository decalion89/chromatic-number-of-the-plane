"""Find where MANY unit circles concur, by counting multiplicities.

Sampling pairs at random found candidates with at most thirteen neighbours,
and growing the graph diluted the high-degree ones rather than adding them.
The sampling was the problem, and the fix falls out of the arithmetic.

A point with k neighbours is the intersection of the unit circles about each
of C(k, 2) PAIRS of them.  So it does not appear once in the candidate list,
it appears C(k,2) times -- a 13-neighbour point appears 78 times, a
20-neighbour point 190.  Enumerate every intersection and count multiplicity,
and the high-concurrence points sort themselves to the top.  Random sampling
throws that signal away.

It is also cheap.  An intersection is constructible only when the pair sits at
a squared distance D with sqrt(D) and sqrt(4-D) both in the field, so only
pairs at RATIONAL squared distance can contribute, and G has about fifty-two
thousand of those rather than the 1.25 million pairs it has in total.  Every
one can be done exactly.

Reported: the multiplicity distribution, the implied neighbour counts, and
then the exact placeability test on the very best ones.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from hn.degrey import build_G, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]


def rsqrt_rational(v):
    """sqrt of a positive rational, as a field element, if it is there."""
    if v <= 0:
        return None
    for r in CLASSES:
        q = v / r
        n2, d2 = q.numerator, q.denominator
        rn, rd = int(round(n2 ** .5)), int(round(d2 ** .5))
        if rn * rn == n2 and rd * rd == d2:
            s = K.rational(Fr(rn, rd))
            return s if r == 1 else K.sqrt(r) * s
    return None


def run(name, P):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    n = len(P)

    # every pair at a rational squared distance, grouped by that distance
    bydist = defaultdict(list)
    for i in range(n - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            v = int(sq[off, 0])
            if v:
                bydist[Fr(v, D2)].append((i, i + 1 + int(off)))
    print(f"{name}: {n} pts, {sum(len(v) for v in bydist.values())} pairs at "
          f"a rational squared distance over {len(bydist)} classes"
          f"  [{time.time()-t0:.0f}s]", flush=True)

    half = K.rational(Fr(1, 2))
    # Points ALREADY in the graph turn up as intersections of their own
    # neighbours -- the hub, degree 60, is the intersection of C(60,2) = 1770
    # pairs -- and they are trivially "placeable" because their own colour is
    # free.  The first run reported a 60-neighbour point for exactly that
    # reason.  Only new points are candidates.
    inP = set(P)
    mult = Counter()
    where = {}
    usable = 0
    for D, pairs in bydist.items():
        if float(D) >= 4 or float(D) <= 0:
            continue
        sD = rsqrt_rational(D)
        s4 = rsqrt_rational(Fr(4) - D)
        if sD is None or s4 is None:
            continue
        usable += 1
        c = s4 * (K.rational(1) / sD) * half
        for i, j in pairs:
            A, B = P[i], P[j]
            mx, my = (A.x + B.x) * half, (A.y + B.y) * half
            nx, ny = -(B.y - A.y) * c, (B.x - A.x) * c
            for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
                if q in inP:
                    continue
                mult[q] += 1
                where[q] = q
    print(f"   {usable} usable classes, {len(mult)} distinct NEW "
          f"intersection points (existing vertices excluded)"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    top = mult.most_common(60)
    print(f"   top multiplicities: {[m for _, m in top[:20]]}", flush=True)

    # exact neighbour counts and the placeability test on the best
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    cls = [[1 + v * k + col for col in range(k)] for v in range(n)]
    for a, b in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + b * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    assert sv.solve()
    cand = [q for q, m in top]
    gb = IntBasis.covering(P + cand)
    gdim, gD2 = gb.dim, gb.D * gb.D
    Prows, crows = gb.rows(P), gb.rows(cand)
    best = 0
    for t, o in enumerate(crows):
        d = Prows - o
        sq = gb._field_square(d[:, :gdim]) + gb._field_square(d[:, gdim:])
        hit = sq[:, 0] == gD2
        for m in range(1, gdim):
            hit &= sq[:, m] == 0
        nb = np.nonzero(hit)[0]
        best = max(best, len(nb))
        if len(nb) >= k and not sv.solve(
                assumptions=[-(1 + int(v) * k) for v in nb]):
            print(f"*** {name}: a point with {len(nb)} neighbours has NO free "
                  f"colour -- chi >= 6 ***", flush=True)
            return True
        if t < 12:
            print(f"      multiplicity {top[t][1]:4d} -> {len(nb)} actual "
                  f"neighbours, placeable", flush=True)
    print(f"   largest neighbourhood among the top 60: {best}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    return False


for nm, P in (("G", build_G(K, as_graph=False)), ("Y", build_Y(K))):
    if run(nm, P):
        break
print("DONE", flush=True)
