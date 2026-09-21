"""The dense closure's middle-sized distance classes.

Two things point the same way.  Sa's three carrying classes sit 5th, 6th and
9th of fourteen by size -- never the biggest.  And G's one hard class, D = 4/9,
is the fifth of its thirty-six and the SMALLEST of the five it was asked.  The
populous classes forbid so much that the solver finds a colouring immediately;
the useful ones are in the middle.

The densest object built here is G closed dihedrally about its own hub: 11047
points at 5.424 edges per vertex, against 4.98 for every bite and thickening.
Its biggest class, D = 1/3 with 45672 pairs, colours instantly, exactly as the
pattern predicts.  Its middle is untested.

Order the closable classes by size and take the middle band, smallest first so
the cheap ones answer before the expensive ones start.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
G = build_G(F, as_graph=False)
C = G[0]
rot = _rot60(F).about(C)
seen, P = set(), []
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, C.y + (C.y - q.y))
            if q not in seen:
                seen.add(q)
                P.append(q)
g = build_graph(P)
E = set((min(a, b), max(a, b)) for a, b in g.edges())
print(f"closure about the hub: {g.n} pts, {len(E)} edges, "
      f"{len(E)/g.n:.3f} per vertex  [{time.time()-t0:.0f}s]", flush=True)

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
clo = sorted((d for d in byd if closable_distance(d)), key=lambda d: len(byd[d]))
print(f"{len(clo)} closable classes, sizes {[len(byd[d]) for d in clo]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

base = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
for a, b in sorted(E):
    for c in range(k):
        base.append([-(1 + a * k + c), -(1 + b * k + c)])
for d in clo:
    pr = byd[d]
    if len(pr) < 40:
        continue
    cls = list(base)
    for i, j in pr:
        for c in range(k):
            cls.append([-(1 + i * k + c), -(1 + j * k + c)])
    t1 = time.time()
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    verdict = ("colours" if ok else
               "*** DOES NOT COLOUR -- THE WEAK PROPERTY HOLDS ***")
    print(f"   D={str(d):9s} {len(pr):7d} pairs -> {verdict} in "
          f"{time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        break
print("DONE", flush=True)
