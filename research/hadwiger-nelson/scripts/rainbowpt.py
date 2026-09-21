"""Ask each candidate point the exact question, instead of counting edges.

Growing by "most new edges" raised density from 4.98 to 5.56 on G without
moving anything, which is the same thing that happened one level down.  Edges
are a proxy and there is no need for a proxy here.

A new point p can be added to a 5-colouring exactly when some colour is free
in its neighbourhood.  So G + p is NOT 5-colourable precisely when, in EVERY
proper 5-colouring of G, the neighbours of p already use all five colours --
when N(p) is rainbow in every colouring.  And that is five SAT calls with
assumptions, one per colour: is there a colouring in which no neighbour of p
takes colour c?  Five UNSATs and the point has nowhere to go.

The calls are assumptions on one bootstrapped solver, so each candidate costs
five incremental solves and nothing is rebuilt.  A candidate needs at least
five neighbours to stand a chance, and the constructible ones here reach
thirteen.

This is exact and it is the question itself, not a correlate of it.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
rng = random.Random(13)
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


def run(name, P, k=5, tries=200000):
    basis0 = IntBasis.covering(P)
    rows0 = basis0.rows(P)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis0, rows0)))
    n = len(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    assert sv.solve(), f"{name} is not {k}-colourable"
    cand = candidates(P, tries)
    print(f"\n{name}: {n} pts, {len(E)} edges, {len(cand)} constructible "
          f"candidates  [{time.time()-t0:.0f}s]", flush=True)
    gb = IntBasis.covering(P + cand)
    dim, D2 = gb.dim, gb.D * gb.D
    Prows = gb.rows(P)
    crows = gb.rows(cand)
    # Colour symmetry makes this ONE call, not five: if some colouring leaves
    # a colour free on N(p), permuting colours leaves colour 0 free, so
    # testing c = 0 decides it.  The five-call version broke at the first free
    # colour anyway, which is why it always reported "0 of 5 blocked" -- the
    # count can only ever be 0 or 5.  One call per candidate means twenty
    # times as many candidates for the same price.
    from collections import Counter
    sizes = Counter()
    placeable = blockedn = 0
    worst = (0, None)
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
        if len(nb) > worst[0]:
            worst = (len(nb), t)
        ass = [-(1 + int(v) * k + c0) for v in nb for c0 in (0,)]
        if sv.solve(assumptions=ass):
            placeable += 1
        else:
            blockedn += 1
            print(f"*** candidate {t}: {len(nb)} neighbours and NO free "
                  f"colour -- {name} + that point is not {k}-colourable, "
                  f"so chi >= {k+1} ***", flush=True)
            return cand[t]
    print(f"   neighbour counts: {dict(sorted(sizes.items()))}", flush=True)
    print(f"   {placeable} candidates with >= {k} neighbours, every one "
          f"placeable; largest neighbourhood {worst[0]}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    return None


for nm, P in (("G", build_G(K, as_graph=False)),
              ("Y", build_Y(K)), ("Sa", build_Sa(K))):
    got = run(nm, P)
    if got is not None:
        break
print("\nDONE", flush=True)
