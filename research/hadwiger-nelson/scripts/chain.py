"""Chain the two: grow G, then ask the exact question of the grown graph.

The exact test says G's largest constructible neighbourhood is thirteen and
every candidate with five or more neighbours is placeable.  It also says what
would change that: bigger neighbourhoods, which need more points, because a
candidate's neighbours beyond the two that construct it are concurrences of
unit circles and concurrences need circles to concur.

The growth is producing exactly that -- 1581 points to 2601 and climbing, 4.98
edges per vertex to 6.15 -- so run the test on the grown graph rather than on
G, at each stage, and watch the largest neighbourhood.  If it climbs past
thirteen the two lines are feeding each other; if it sticks, the ceiling is
structural and worth knowing as one.

Every candidate that fails is a proper five-colouring found, which is a real
answer.  Any candidate that succeeds is a point with no colour available, and
the graph plus that point needs six.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
rng = random.Random(29)
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


def candidates(P, tries):
    half = K.rational(Fr(1, 2))
    out, seen = [], set(P)
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
            if q not in seen:
                seen.add(q)
                out.append(q)
    return out


def graph(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    return b, r, sorted(set((min(a, c), max(a, c))
                            for a, c in fast_edges_complete(b, r)))


def probe(P, tag, tries=120000):
    b, r, E = graph(P)
    n = len(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {tag}: {n} pts NOT 5-COLOURABLE", flush=True)
        return True
    cand = candidates(P, tries)
    gb = IntBasis.covering(P + cand)
    dim, D2 = gb.dim, gb.D * gb.D
    Prows, crows = gb.rows(P), gb.rows(cand)
    sizes, worst, placeable = Counter(), 0, 0
    for t, o in enumerate(crows):
        d = Prows - o
        sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
        hit = sq[:, 0] == D2
        for m in range(1, dim):
            hit &= sq[:, m] == 0
        nb = np.nonzero(hit)[0]
        sizes[len(nb)] += 1
        if len(nb) < k:
            continue
        worst = max(worst, len(nb))
        if sv.solve(assumptions=[-(1 + int(v) * k) for v in nb]):
            placeable += 1
        else:
            print(f"*** {tag}: a point with {len(nb)} neighbours has NO free "
                  f"colour -- chi >= 6 ***  [{time.time()-t0:.0f}s]",
                  flush=True)
            return True
    sv.delete()
    big = {a: c for a, c in sorted(sizes.items()) if a >= 5}
    print(f"  {tag}: {n} pts, {len(E)} edges, {len(E)/n:.2f} per vertex, "
          f"{len(cand)} candidates; >=5 neighbours {big}; largest {worst}, "
          f"all {placeable} placeable  [{time.time()-t0:.0f}s]", flush=True)
    return False


P = build_G(K, as_graph=False)
have = set(P)
probe(P, "G")
BATCH = 120
for step in range(1, 12):
    cand = [c for c in candidates(P, 60000) if c not in have]
    if not cand:
        break
    gb = IntBasis.covering(P + cand[:1500])
    dim, D2 = gb.dim, gb.D * gb.D
    Prows = gb.rows(P)
    scored = []
    for c in cand[:1500]:
        o = gb.rows([c])[0]
        d = Prows - o
        sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
        hit = sq[:, 0] == D2
        for m in range(1, dim):
            hit &= sq[:, m] == 0
        scored.append((int(hit.sum()), c))
    scored.sort(key=lambda t: -t[0])
    add = [c for s, c in scored[:BATCH] if s >= 3]
    if not add:
        break
    P = list(dict.fromkeys(P + add))
    have = set(P)
    if probe(P, f"grown x{step}"):
        break
print("DONE", flush=True)
