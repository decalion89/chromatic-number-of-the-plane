"""Drop the antipodal requirement: all pairs at one closable distance.

Antipodality was de Grey's SYMMETRIC choice, not the lemma's requirement.  The
weak property needs a set of non-adjacent pairs, all at the same distance d,
such that forbidding every one of them kills colourability -- then in every
proper colouring at least one of them is monochromatic, and they all share the
distance, so they all share the spindle rotation ANGLE.  Only the centre of
that rotation differs from pair to pair.

That is a far larger family than the antipodal one.  G has 21344 non-edge
pairs at a rational closable distance, spread over 36 distances, against 119
antipodal pairs across 8 rings about its hub.  The D = 1/3 class alone holds
6510 pairs, and forbidding six thousand pairs is a serious constraint where
forbidding three is not.

Thirty-six calls, one per closable distance, and the same one-call identity:
forbidding a pair from being monochromatic is adding the edge.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def by_distance(P):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    out = defaultdict(list)
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E:
                out[Fr(int(sq[off, 0]), D2)].append((i, j))
    return g.n, E, out


def scan(name, P, k):
    n, E, byd = by_distance(P)
    base = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            base.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv0 = Solver(name="cd15", bootstrap_with=base)
    assert sv0.solve(), f"{name} is not {k}-colourable"
    sv0.delete()
    clo = sorted((d for d in byd if closable_distance(d)),
                 key=lambda d: -len(byd[d]))
    print(f"\n{name} at {k}: {n} pts, {len(E)} edges, {len(byd)} rational "
          f"distance classes, {len(clo)} closable  [{time.time()-t0:.0f}s]",
          flush=True)
    hits = []
    for d in clo:
        pr = byd[d]
        cls = list(base)
        for i, j in pr:
            for c in range(k):
                cls.append([-(1 + i * k + c), -(1 + j * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        mark = "" if ok else "  *** UNCOLOURABLE: SOME PAIR ALWAYS MONO ***"
        if not ok:
            hits.append((d, len(pr)))
        print(f"   D={str(d):8s} {len(pr):6d} pairs forbidden -> "
              f"{'colours' if ok else 'DOES NOT COLOUR'}{mark}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"  {len(hits)} of {len(clo)} closable distance classes carry it",
          flush=True)
    return hits


scan("Sa", build_Sa(F), 4)
scan("G", build_G(F, as_graph=False), 5)
print("\nDONE", flush=True)
