"""The by-distance gateway on the biggest five-chromatic objects here.

Sa carries the weak property on three closable distance classes and the
carrying ones are not the populous ones -- 4/9, 4, 16/9 against 1/3, 7/3, 3,
5/9 which all colour.  The distances themselves are 2/3, 2 and 4/3: an
arithmetic progression of step 2/3, which is Sa's own scale.

G fails every antipodal ring at five and the by-distance scan is running on
it.  What has not been asked is whether anything BIGGER in the same shape
carries it, so ask the two most Sa-like objects available.

  G closed about its hub: 11047 points, six-fold symmetric about a degree-60
  vertex, which is exactly the shape Sa has.
  G bitten at the sqrt(17) ring: the union with the most coupling built here,
  48 ring points against de Grey's 6.

Top classes only -- the call cost grows with the number of pairs forbidden and
the populous classes are the ones Sa says do not matter anyway, so take a
spread rather than the top of the list.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.field import Field, embed
from hn.geometry import (DEGREY_FIELD as F, Point, _rot60, rotation_joining,
                         Rotation)
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()


def scan(name, P, limit=18):
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
    sv0 = Solver(name="cd15", bootstrap_with=base)
    ok = sv0.solve()
    sv0.delete()
    clo = sorted((d for d in byd if closable_distance(d)),
                 key=lambda d: -len(byd[d]))
    print(f"\n{name}: {g.n} pts, {len(E)} edges, "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}, "
          f"{len(clo)} closable classes  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        return
    # a spread: the biggest few and the smallest few, since Sa says size is
    # not what decides
    take = clo[:limit // 2] + clo[-(limit // 2):]
    for d in dict.fromkeys(take):
        pr = byd[d]
        cls = list(base)
        for i, j in pr:
            for c in range(k):
                cls.append([-(1 + i * k + c), -(1 + j * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        good = sv.solve()
        sv.delete()
        print(f"   D={str(d):9s} {len(pr):7d} pairs -> "
              f"{'colours' if good else '*** DOES NOT COLOUR ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)


G = build_G(F, as_graph=False)
C = G[0]
rot = _rot60(F).about(C)
seen, H = set(), []
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
                H.append(q)
scan("G closed about its hub", H)

K17 = Field((3, 5, 7, 11, 17))
GB = [Point(embed(p.x, K17), embed(p.y, K17)) for p in G]
CB = GB[0]
b = rotation_joining(Fr(5, 3), K17)
f, iv = b.about(CB), Rotation(b.cos, -b.sin).about(CB)
seen, U = set(), []
for p in list(GB) + [f(p) for p in GB] + [iv(p) for p in GB]:
    if p not in seen:
        seen.add(p)
        U.append(p)
scan("G bitten at the sqrt(17) ring, both ways", U)
print("\nDONE", flush=True)
