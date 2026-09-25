"""Combine the carrying classes, and ask five colours again.

Four classes carry the weak property for Sa at four colours -- 4/9 with 393
pairs, 4 with 273, 16/9 with 156, and 16 with just three, which are de Grey's
antipodal pairs.  At five colours none of them carries alone, but 4/9 and 16/9
are the two that behave unlike all the others: 10s for Sa, 577s for Y, 372s
for the hub closure, over an hour for G, against none for every other class.

Forbidding two classes at once is a weaker statement than forbidding one --
"some pair at distance 2/3 OR 4/3 is monochromatic in every 5-colouring"
rather than "some pair at 2/3 is".  It is still a real property, and a
multispindle over two distances is more awkward than over one rather than
impossible.  And it is one SAT call.

Ordered small graphs first, and subsets before the full set, so the cheap
answers arrive before the expensive ones start.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.fast import IntBasis
from pysat.solvers import Solver

t0 = time.time()
CARRY = [Fr(4, 9), Fr(16, 9), Fr(4), Fr(16)]


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
    return g.n, E, byd


def ask(name, n, E, byd, ds, k):
    pr = [p for d in ds for p in byd.get(d, [])]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E) + pr:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    t1 = time.time()
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    tag = " + ".join(str(d) for d in ds)
    print(f"  {name} at {k}, {{{tag}}}: {len(pr)} pairs -> "
          f"{'colours' if ok else '*** DOES NOT COLOUR -- CARRIES IT ***'}"
          f" in {time.time()-t1:.0f}s  [{time.time()-t0:.0f}s]", flush=True)
    return not ok


for name, P in (("Sa", build_Sa(F)), ("Y", build_Y(F))):
    n, E, byd = classes(P)
    print(f"\n{name}: {n} pts  [{time.time()-t0:.0f}s]", flush=True)
    for r in (2, 3, 4):
        for ds in itertools.combinations(CARRY, r):
            if ask(name, n, E, byd, ds, 5):
                print(f"  *** {name} carries it at FIVE on {ds} ***",
                      flush=True)
                break
print("\nDONE", flush=True)
