"""Grow G itself, which is already five-chromatic, toward six.

The seed programme attacked the wrong level.  Growing from Sa -- four
chromatic -- produced a better SEED, five carrying classes instead of four,
and nothing at five: the levels decouple, so a better seed does not climb.

The symmetric move was never made.  G is 5-chromatic already.  Grow IT the
same way: add points at distance one from two points already present, which
is where two unit circles meet and is constructible when sqrt(D) and
sqrt(4-D) both lie in the field, choosing each batch to create the most new
edges, and ask after every batch whether the graph still colours with five.

Two differences from the seed version, both forced by what G is.  G is not
D6-symmetric about the origin -- its symmetry is about the spindle pivot
(-2, 0) -- so points are added singly rather than in orbits.  And single
points move a 1581-vertex graph slowly, so they go in batches.

This is the direct attack: no pipeline, no seed, no spindle to follow.  Just
whether a five-chromatic unit-distance graph can be pushed off five by adding
the points the field allows.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
rng = random.Random(11)
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


def candidates(P, tries=40000):
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
            if q not in seen:
                seen.add(q)
                out.append(q)
    return out


def graph(P):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    return basis, rows, E


def kcol(n, E, kk):
    cls = [[1 + v * kk + c for c in range(kk)] for v in range(n)]
    for a, b in E:
        for c in range(kk):
            cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


P = build_G(K, as_graph=False)
have = set(P)
basis, rows, E = graph(P)
print(f"G: {len(P)} pts, {len(E)} edges, {len(E)/len(P):.2f} per vertex, "
      f"5-colourable {kcol(len(P), E, 5)}  [{time.time()-t0:.0f}s]",
      flush=True)

BATCH = 60
for step in range(1, 26):
    cand = [c for c in candidates(P) if c not in have]
    if not cand:
        print("   no constructible candidates", flush=True)
        break
    # rank candidates by how many unit edges each would bring, vectorised
    gb = IntBasis.covering(P + cand[:1200])
    dim, D2 = gb.dim, gb.D * gb.D
    Prows = gb.rows(P)
    scored = []
    for c in cand[:1200]:
        o = gb.rows([c])[0]
        d = Prows - o
        sq = gb._field_square(d[:, :dim]) + gb._field_square(d[:, dim:])
        hit = sq[:, 0] == D2
        for m in range(1, dim):
            hit &= sq[:, m] == 0
        scored.append((int(hit.sum()), c))
    scored.sort(key=lambda t: -t[0])
    add = [c for s, c in scored[:BATCH] if s >= 2]
    if not add:
        print("   nothing worth adding", flush=True)
        break
    P = list(dict.fromkeys(P + add))
    have = set(P)
    basis, rows, E = graph(P)
    ok = kcol(len(P), E, 5)
    print(f"  step {step}: +{len(add)} -> {len(P)} pts, {len(E)} edges, "
          f"{len(E)/len(P):.2f} per vertex, "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  (best candidate brought {scored[0][0]} edges)"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
