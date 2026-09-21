"""The whole family at five colours, class by class.

The conclusion the measurements converge on is that the bite sharpens a weak
property rather than creating one, and that G carries none to sharpen.  What
closes that statement is the rest of the family: Sa and Y are the two graphs
that DO carry it at four, so whether they carry anything at five is the
cleanest test of whether the property is about the graph or about the number
of colours.

They are small -- 397 and 791 points -- so every closable class is affordable,
and the answer takes minutes rather than the hours G's classes cost.

If both come back empty, the family is uniformly empty at five: the two graphs
that carry the property at four colours lose it at five, and the graph built
from them by spindling never had it.  That is a closing statement rather than
another negative.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y, build_Sb
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def scan(name, P, k):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    byd = defaultdict(list)
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E:
                byd[Fr(int(sq[off, 0]), D2)].append((i, j))
    base = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E):
        for c in range(k):
            base.append([-(1 + a * k + c), -(1 + b * k + c)])
    clo = sorted((d for d in byd if closable_distance(d)),
                 key=lambda d: len(byd[d]))
    hits = []
    for d in clo:
        pr = byd[d]
        cls = list(base)
        for i, j in pr:
            for c in range(k):
                cls.append([-(1 + i * k + c), -(1 + j * k + c)])
        t1 = time.time()
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        if not ok:
            hits.append((d, len(pr)))
            print(f"   *** {name} at {k}, D={d}: {len(pr)} pairs, DOES NOT "
                  f"COLOUR ***  [{time.time()-t0:.0f}s]", flush=True)
        elif time.time() - t1 > 5:
            print(f"   {name} at {k}, D={d}: {len(pr)} pairs, colours but "
                  f"took {time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]",
                  flush=True)
    print(f"{name} at {k}: {g.n} pts, {len(clo)} closable classes, "
          f"{len(hits)} carry the weak property"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    return hits


scan("Sa", build_Sa(F), 5)
scan("Sb", build_Sb(F), 5)
scan("Y", build_Y(F), 5)
print("\n-- and the four-colour controls, which must come back non-empty --",
      flush=True)
scan("Sa", build_Sa(F), 4)
print("\nDONE", flush=True)
