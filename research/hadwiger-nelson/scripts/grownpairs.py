"""The exact forced-pair census on the GROWN graphs, which never got it.

The bipartite theorem says a five-colour obstruction cannot be local: every
neighbourhood is 2-chromatic on its own and, measured, every one of G's 1581
vertices and 400 best new candidates is forced to exactly two colours.  So the
only mechanism left is de Grey's, which is non-local by construction -- a PAIR
monochromatic in every colouring, spindled about one end.

That census was run on G, Y, Sa and their closures and bites, and came back
empty at five: 21344 closable non-edge pairs of G, 0 forced, complete in
twenty-one seconds.  It has never been run on the graphs the growth produced,
which are the densest objects here -- G grown to three thousand points at 6.35
edges per vertex against its own 4.98.

Same test, same exactness: every non-edge pair at a rational squared distance
the field can spindle, one SAT call each, forced-same iff the graph plus that
edge stops colouring.  Growing the graph can only add forced pairs, never
remove them, since forcing is monotone under adding vertices -- so if the
growth has bought anything at all, this is where it shows.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, forced_same
from pysat.solvers import Solver

k = 5
t0 = time.time()
rng = random.Random(41)
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


def candidates(P, tries=60000):
    half = K.rational(Fr(1, 2))
    out, seen = [], set(P)
    n = len(P)
    for _ in range(tries):
        A, B = P[rng.randrange(n)], P[rng.randrange(n)]
        if A == B:
            continue
        D = A.dist2(B)
        f = float(D)
        if not .05 < f < 3.99:
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


def census(name, P):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    n = len(P)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    Eset = set(E)
    pairs, byd = [], Counter()
    for i in range(n - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) in Eset:
                continue
            dd = Fr(int(sq[off, 0]), D2)
            if dd and closable_distance(dd):
                pairs.append((i, j))
                byd[dd] += 1
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"*** {name}: {n} pts NOT 5-COLOURABLE ***", flush=True)
        return True
    forced = [(i, j) for i, j in pairs if forced_same(sv, i, j, k)]
    sv.delete()
    print(f"  {name}: {n} pts, {len(E)} edges, {len(E)/n:.2f} per vertex, "
          f"{len(pairs)} closable non-edge pairs over {len(byd)} classes, "
          f"{len(forced)} FORCED  [{time.time()-t0:.0f}s]", flush=True)
    if forced:
        print(f"  *** forced pairs: {forced[:5]} -- spindle them ***",
              flush=True)
        return True
    return False


P = build_G(K, as_graph=False)
have = set(P)
census("G", P)
BATCH = 150
for step in range(1, 11):
    cand = [c for c in candidates(P) if c not in have]
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
    if census(f"grown x{step}", P):
        break
print("DONE", flush=True)
