"""Grow from Sa, which already works, instead of from nothing.

Growing a D6-symmetric set from the unit orbit stalls at 31 points and chi = 3:
the candidates are unit-circle intersections, those need sqrt(D) and sqrt(4-D)
both in the field, and after three steps every pair is at an irrational
squared distance the test cannot handle.  Widening the field from four
radicals to six changed nothing, which says the obstruction is the test and
the starting point, not the arithmetic.

Sa is the starting point that works.  It is 4-chromatic, it carries the weak
property on four closable classes, and its pairwise distances include plenty
of rational ones -- 1/3 with 38 pairs, 1/9 with 21, 1 with 18 -- so
constructible candidates exist.  Growing from there keeps everything Sa has
and asks what one more orbit buys.

Each step: find the points at distance one from two points of Sa, take the
orbit of whichever adds the most edges, and measure chi and the weak property
at five.  If chi ever reaches five while the D6 symmetry holds, that is a
five-chromatic rotationally symmetric graph -- the object the whole structural
account says is needed, and one the family does not contain.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()
W = _rot60(K)
rng = random.Random(7)
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]


def rsqrt(e):
    if any(x for x in e.c[1:]):
        return None
    v = e.c[0]
    if v <= 0:
        return None
    for r in CLASSES:
        q = v / r
        n2, d2 = q.numerator, q.denominator
        rn, rd = int(round(n2 ** 0.5)), int(round(d2 ** 0.5))
        if rn * rn == n2 and rd * rd == d2:
            s = K.rational(Fr(rn, rd))
            return s if r == 1 else K.sqrt(r) * s
    return None


def orbit(p):
    out, q = [], p
    for _ in range(6):
        out.append(q)
        out.append(Point(q.x, -q.y))
        q = W(q)
    return out


def edges(P):
    g = build_graph(P)
    return set((min(a, b), max(a, b)) for a, b in g.edges())


def chrom(P, E, hi=5):
    n = len(P)
    for k in range(4, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        if ok:
            return k
    return f">{hi}"


def candidates(P, tries=30000):
    half = K.rational(Fr(1, 2))
    out, seen = [], set()
    n = len(P)
    for _ in range(tries):
        A, B = P[rng.randrange(n)], P[rng.randrange(n)]
        if A == B:
            continue
        D = A.dist2(B)
        f = float(D)
        if f > 3.99 or f < 0.05:
            continue
        sD, s4 = rsqrt(D), rsqrt(K.rational(4) - D)
        if sD is None or s4 is None:
            continue
        inv = K.rational(1) / sD
        mx, my = (A.x + B.x) * half, (A.y + B.y) * half
        nx = -(B.y - A.y) * inv * s4 * half
        ny = (B.x - A.x) * inv * s4 * half
        for p in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if p not in seen and float(p.norm2()) < 9:
                seen.add(p)
                out.append(p)
    return out


P = build_Sa(K)
E = edges(P)
print(f"Sa: {len(P)} pts, {len(E)} edges, chi = {chrom(P, sorted(E))}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
have = set(P)
for step in range(1, 40):
    cand = [c for c in candidates(P) if c not in have]
    if not cand:
        print("   no constructible candidates", flush=True)
        break
    # Gain, vectorised.  The first version compared exact Points in a
    # Python loop -- 250 candidates x 12 orbit points x 400 set points is
    # over a million Fraction comparisons, and it was SLOWER than rebuilding
    # the graph, which at least used the int64 edge finder.  int64 rows and
    # one field-square per orbit point is the right shape.
    pool = P + [q for c in cand[:250] for q in orbit(c)]
    gb = IntBasis.covering(pool)
    dim, D2 = gb.dim, gb.D * gb.D
    Prows = gb.rows(P)

    def gain_of(c):
        orb = [q for q in dict.fromkeys(orbit(c)) if q not in have]
        if not orb:
            return -1
        orows = gb.rows(orb)
        g2 = 0
        for o in orows:
            d = Prows - o
            sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
            hit = sq[:, 0] == D2
            for m in range(1, dim):
                hit &= sq[:, m] == 0
            g2 += int(hit.sum())
        for a in range(len(orows)):
            d = orows[a + 1:] - orows[a]
            if not len(d):
                continue
            sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
            hit = sq[:, 0] == D2
            for m in range(1, dim):
                hit &= sq[:, m] == 0
            g2 += int(hit.sum())
        return g2

    best, gain = None, -1
    for c in cand[:250]:
        g2 = gain_of(c)
        if g2 > gain:
            best, gain = c, g2
    P = list(dict.fromkeys(P + orbit(best)))
    have = set(P)
    E = edges(P)
    k = chrom(P, sorted(E))
    print(f"  step {step}: {len(P)} pts, {len(E)} edges, "
          f"{len(E)/len(P):.2f} per vertex, chi = {k}  (gain {gain}, "
          f"{len(cand)} candidates)  [{time.time()-t0:.0f}s]", flush=True)
    if k == ">5":
        print("*** chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
