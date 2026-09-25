"""Not "is it blocked" but "how many colours does its neighbourhood force".

The placeability test is a yes/no and every candidate answers no.  The same
machinery gives a graded and still EXACT answer: for a candidate p, what is
the minimum, over all proper 5-colourings of G, of the number of distinct
colours appearing on N(p)?

    minimum 5  ->  N(p) is rainbow in every colouring, p has no colour, chi>=6
    minimum 4  ->  one colour away
    minimum 2  ->  nowhere near

It is computed by assumptions and nothing else.  "N(p) uses at most j colours"
is, by colour symmetry, "no vertex of N(p) takes any of the colours
j, j+1, ..., 4" -- so one solve per j and the smallest feasible j is the
answer.  Four calls per candidate, on a solver bootstrapped once.

That turns a wall of identical negatives into a distribution, and a
distribution says whether the ceiling of thirteen neighbours is the binding
constraint or a red herring: a five-point neighbourhood forced rainbow would
win just as well as a thirteen-point one, and the minimum is what decides.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]


def rsqrt_rational(v):
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


def run(name, P, keep=400):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    n = len(P)
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
    half = K.rational(Fr(1, 2))
    inP = set(P)
    mult = Counter()
    for D, pairs in bydist.items():
        if not 0 < float(D) < 4:
            continue
        sD, s4 = rsqrt_rational(D), rsqrt_rational(Fr(4) - D)
        if sD is None or s4 is None:
            continue
        c = s4 * (K.rational(1) / sD) * half
        for i, j in pairs:
            A, B = P[i], P[j]
            mx, my = (A.x + B.x) * half, (A.y + B.y) * half
            nx, ny = -(B.y - A.y) * c, (B.x - A.x) * c
            for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
                if q not in inP:
                    mult[q] += 1
    cand = [q for q, _ in mult.most_common(keep)]
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    cls = [[1 + v * k + col for col in range(k)] for v in range(n)]
    for a, b in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + b * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    assert sv.solve()
    gb = IntBasis.covering(P + cand)
    gdim, gD2 = gb.dim, gb.D * gb.D
    Prows, crows = gb.rows(P), gb.rows(cand)
    dist = Counter()
    bysize = defaultdict(Counter)
    for t, o in enumerate(crows):
        d = Prows - o
        sq = gb._field_square(d[:, :gdim]) + gb._field_square(d[:, gdim:])
        hit = sq[:, 0] == gD2
        for m in range(1, gdim):
            hit &= sq[:, m] == 0
        nb = [int(v) for v in np.nonzero(hit)[0]]
        if len(nb) < 2:
            continue
        lo = k
        for j in range(1, k + 1):
            ass = [-(1 + v * k + col) for v in nb for col in range(j, k)]
            if sv.solve(assumptions=ass):
                lo = j
                break
        dist[lo] += 1
        bysize[len(nb)][lo] += 1
        if lo >= k:
            print(f"*** {name}: candidate {t} with {len(nb)} neighbours "
                  f"forces all {k} colours -- chi >= {k+1} ***", flush=True)
            return True
    print(f"\n{name}: {n} pts, {len(cand)} top candidates tested", flush=True)
    print(f"   minimum colours forced on N(p): {dict(sorted(dist.items()))}",
          flush=True)
    for size in sorted(bysize)[-6:]:
        print(f"      {size:3d} neighbours -> {dict(sorted(bysize[size].items()))}",
              flush=True)
    print(f"   [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    return False


for nm, P in (("G", build_G(K, as_graph=False)), ("Y", build_Y(K)),
              ("Sa", build_Sa(K))):
    if run(nm, P):
        break
print("DONE", flush=True)
