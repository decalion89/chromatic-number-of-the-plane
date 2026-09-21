"""Settle the grown seed at five colours with the calibrated local search.

CDCL has been on the first of Sa''s closable classes for five minutes.  That
is not evidence: twice in this session a hard instance turned out satisfiable
after hours of it, on G's class 4/9 and on Sa with all four carrying classes,
and the rule that came out of those is that only LANDING proves anything.

TabuCol, calibrated in both directions earlier -- it reaches zero in seconds
on instances known satisfiable and plateaus at 37 on one known unsatisfiable
-- settled both of those in minutes.  So put it on every closable class of the
grown seed at five colours, smallest first, and let CDCL keep the ones it does
not land.

A landing rules a class out, fast.  A plateau rules nothing in.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from tabucol import tabucol
import random

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
            sq = K.rational(Fr(rn, rd))
            return sq if r == 1 else K.sqrt(r) * sq
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
        for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if q not in seen and float(q.norm2()) < 9:
                seen.add(q)
                out.append(q)
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
print(f"grown seed: {len(P)} pts, {len(E)} edges  [{time.time()-t0:.0f}s]",
      flush=True)

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
clo = sorted((d for d in byd if closable_distance(d)),
             key=lambda d: len(byd[d]))
print(f"{len(clo)} closable classes, sizes "
      f"{[len(byd[d]) for d in clo[:14]]}  [{time.time()-t0:.0f}s]",
      flush=True)

for d in clo[:14]:
    arr = np.array(E + byd[d], dtype=np.int64)
    got = []
    for sd in range(4):
        b, it = tabucol(len(P), arr, 5, seed=300 + sd, iters=600_000)
        got.append(b)
        if b == 0:
            break
    lo = min(got)
    print(f"   D={str(d):8s} {len(byd[d]):5d} pairs -> "
          f"{'COLOURS (landed)' if lo == 0 else f'no landing, best {lo}'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
