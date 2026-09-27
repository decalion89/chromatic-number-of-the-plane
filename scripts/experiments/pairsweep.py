"""Every PAIR of closable classes, which has never been swept.

Single classes have been done -- Sa carries four of fourteen at four colours
and none at five -- and so has the union of all four carrying ones, which
colours.  What has not been asked is the ninety-one PAIRS.

A pair is a real statement: "in every 5-colouring some pair of points is
monochromatic at distance d1 or d2".  It is weaker than a single class and
stronger than the four-class union, and it sits in the gap where nothing has
been measured.  Ninety-one SAT calls on 397 points, so the whole sweep is
affordable where the four-class instance alone took 850 seconds.

Ordered by combined size, smallest first, since the four-class run showed the
big unions colour easily and Sa's carrying classes at four were the small and
middle ones, never the biggest.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()


def classes(P):
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
    return g.n, sorted(E), byd


name = sys.argv[1] if len(sys.argv) > 1 else "Sa"
P = build_Sa(F) if name == "Sa" else build_Y(F)
n, E, byd = classes(P)
clo = sorted(d for d in byd if closable_distance(d))
print(f"{name}: {n} pts, {len(E)} unit edges, {len(clo)} closable classes "
      f"-> {len(clo)*(len(clo)-1)//2} pairs  [{time.time()-t0:.0f}s]",
      flush=True)
base = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in E:
    for c in range(k):
        base.append([-(1 + a * k + c), -(1 + b * k + c)])
pairs = sorted(itertools.combinations(clo, 2),
               key=lambda t: len(byd[t[0]]) + len(byd[t[1]]))
hits, slow = [], []
for d1, d2 in pairs:
    pr = byd[d1] + byd[d2]
    cls = list(base)
    for i, j in pr:
        for c in range(k):
            cls.append([-(1 + i * k + c), -(1 + j * k + c)])
    t1 = time.time()
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    dt = time.time() - t1
    if not ok:
        hits.append((d1, d2, len(pr)))
        print(f"  *** {d1} + {d2}: {len(pr)} pairs, DOES NOT COLOUR ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    elif dt > 20:
        slow.append((d1, d2, len(pr), round(dt)))
        print(f"  {d1} + {d2}: {len(pr)} pairs, colours but took {dt:.0f}s"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{name}: {len(hits)} of {len(pairs)} pairs carry the weak property "
      f"at five; {len(slow)} took over twenty seconds: {slow}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
