"""Cliques in the combined graph: a K6 would settle it in six vertices.

Forbidding a distance class IS adding those pairs as edges, so the question
"does Sa colour with classes 4/9, 16/9, 4 and 16 forbidden" is just "is the
combined graph 5-colourable".  That graph has a clique number, and a clique
number of six answers it outright with a six-point witness -- no SAT call, no
hours.

Worth asking before waiting: the combined graph has average degree 14, which
is comfortably enough room for a K6 to exist.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis

CARRY = {Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)}
t0 = time.time()


def combined(P, classes):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    extra = []
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E and Fr(int(sq[off, 0]), D2) in classes:
                extra.append((i, j))
    n = g.n
    A = np.zeros((n, n), dtype=bool)
    for a, b in list(E) + extra:
        A[a, b] = A[b, a] = True
    return n, A, len(E), len(extra)


def max_clique(A, n, budget=120.0):
    best = [0]
    t1 = time.time()

    def ext(size, Pset):
        if time.time() - t1 > budget:
            return
        if not Pset.size:
            if size > best[0]:
                best[0] = size
            return
        if size + Pset.size <= best[0]:
            return
        while Pset.size:
            if size + Pset.size <= best[0]:
                return
            v = int(Pset[0])
            ext(size + 1, Pset[A[v, Pset]])
            Pset = Pset[1:]

    order = np.argsort(-A.sum(axis=1)).astype(np.int64)
    ext(0, order)
    return best[0]


for name, P in (("Sa", build_Sa(F)), ("Y", build_Y(F))):
    n, A, ne, nx = combined(P, CARRY)
    deg = A.sum(axis=1)
    w = max_clique(A, n)
    verdict = ("*** K6 -- the combined graph needs six colours ***"
               if w >= 6 else f"lower bound {w}, not enough on its own")
    print(f"{name} + the four classes: {n} pts, {ne} unit edges + {nx} added, "
          f"mean degree {deg.mean():.1f}, max {deg.max()}", flush=True)
    print(f"   largest clique found: {w} -> {verdict}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
