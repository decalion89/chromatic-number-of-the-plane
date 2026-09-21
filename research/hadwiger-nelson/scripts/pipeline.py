"""The whole of de Grey's pipeline, on a grown seed instead of his.

Greedy growth from Sa raises density hard -- 4.97 to 6.24 edges per vertex
over eighteen orbits, past everything in the family -- and never moves chi.
That is not a failure of the growth: Sa is 4-chromatic, and what takes it to
five is not bulk but de Grey's two moves, the bite and the spindle.  So put
the grown seed through them.

  grow Sa by N orbits, keeping D6 symmetry         -> Sa'
  check Sa' still carries the weak property at 4   (Sa carries it on 4/9,
                                                    16/9, 4 and 16)
  bite at the ring that carries it, sharpening     -> Y'
  spindle about one end of the forced pair         -> G'
  ask whether G' carries anything at five

Sa' is denser than Sa and has the same symmetry, so if the property survives
the growth it survives on a better object, and the whole chain shifts.  If the
property does NOT survive, that is the sharper finding: it would say growth
destroys exactly what makes the seed work, which is worth knowing and is not
obvious either way.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, _rot60, rotation_joining
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
        rn, rd = int(round(n2 ** .5)), int(round(d2 ** .5))
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
    return list(dict.fromkeys(out))


def edges(P):
    g = build_graph(P)
    return sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))


def kcol(n, E, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


def classes_of(P, E):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    Eset = set(E)
    byd = defaultdict(list)
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in Eset:
                byd[Fr(int(sq[off, 0]), D2)].append((i, j))
    return byd


def weak(P, E, k, cap=10):
    byd = classes_of(P, E)
    clo = sorted((d for d in byd if closable_distance(d)),
                 key=lambda d: len(byd[d]))
    hits = []
    for d in clo[:cap]:
        if not kcol(len(P), E + byd[d], k):
            hits.append((d, len(byd[d])))
    return len(clo), hits


def candidates(P, tries=20000):
    half = K.rational(Fr(1, 2))
    out, seen = [], set()
    n = len(P)
    for _ in range(tries):
        A, B = P[rng.randrange(n)], P[rng.randrange(n)]
        if A == B:
            continue
        D = A.dist2(B)
        f = float(D)
        if f > 3.99 or f < .05:
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
have = set(P)
for step in range(1, 13):
    cand = [c for c in candidates(P) if c not in have]
    if not cand:
        break
    pool = P + [q for c in cand[:200] for q in orbit(c)]
    gb = IntBasis.covering(pool)
    dim, D2 = gb.dim, gb.D * gb.D
    Prows = gb.rows(P)
    best, gain = None, -1
    for c in cand[:200]:
        orb = [q for q in orbit(c) if q not in have]
        if not orb:
            continue
        g2 = 0
        for o in gb.rows(orb):
            d = Prows - o
            sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
            hit = sq[:, 0] == D2
            for m in range(1, dim):
                hit &= sq[:, m] == 0
            g2 += int(hit.sum())
        if g2 > gain:
            best, gain = c, g2
    P = list(dict.fromkeys(P + orbit(best)))
    have = set(P)
E = edges(P)
print(f"grown seed Sa': {len(P)} pts, {len(E)} edges, "
      f"{len(E)/len(P):.2f} per vertex  [{time.time()-t0:.0f}s]", flush=True)
print(f"   4-colourable: {kcol(len(P), E, 4)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
nclo, hits = weak(P, E, 4)
print(f"   weak property at FOUR: {len(hits)} of the {min(nclo,10)} smallest "
      f"closable classes: {hits}  [{time.time()-t0:.0f}s]", flush=True)
# The grown seed carries the property on FIVE classes where Sa carries four,
# and de Grey's own -- D = 16 with exactly three pairs -- survives untouched.
# The new one is D = 20/3 with 102 pairs.
#
# The first version then applied rotation_joining(D) for the carrying class,
# which is the SPINDLE rotation for that distance, not the bite.  De Grey's
# order is: the class D = 16 is carried, the BITE is rho_4 -- the rotation
# that makes the ring of squared radius 4 touch its own image -- and only
# after it has sharpened "one of three" into a named pair does the spindle at
# rho_16 apply.  Biting with the spindle rotation skips the sharpening.
print(f"\n   weak property at FIVE on the grown seed itself:", flush=True)
n5, h5 = weak(P, E, 5)
print(f"      {len(h5)} of the {min(n5,10)} smallest closable classes: {h5}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

if hits:
    # the bite that sharpens: rho_4 about the origin, de Grey's own
    bite = rotation_joining(Fr(4), K)
    Yp = list(dict.fromkeys(P + [bite(p) for p in P]))
    EY = edges(Yp)
    ok4 = kcol(len(Yp), EY, 4)
    print(f"\n   bitten with rho_4 (the sharpening rotation): {len(Yp)} pts, "
          f"{len(EY)} edges, 4-colourable {ok4}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if ok4:
        # does it now force a NAMED pair at squared distance 16?
        byd = classes_of(Yp, EY)
        pr = byd.get(Fr(16), [])
        cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(Yp))]
        for a, b in EY:
            for c in range(4):
                cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        assert sv.solve()
        forced = [(i, j) for i, j in pr
                  if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)])]
        sv.delete()
        print(f"      {len(pr)} pairs at squared distance 16, {len(forced)} "
              f"FORCED SAME at four colours  [{time.time()-t0:.0f}s]",
              flush=True)
        if forced:
            print(f"      *** a named forced pair: {forced[:3]} -- the "
                  f"spindle applies ***", flush=True)
        n5b, h5b = weak(Yp, EY, 5)
        print(f"      weak property at FIVE on the bitten graph: {len(h5b)} "
              f"of {min(n5b,10)}: {h5b}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
