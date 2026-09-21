"""How many closable classes must be unioned before G stops colouring?

Single classes were tested and every one colours; the four carrying classes of
Sa were tested together and they colour too.  What was never asked of G is the
graded version: forbid the classes one at a time, cumulatively, and count how
many it takes.

That number is a real measure and it is exact.  If the union of two or three
classes stops colouring, then in every 5-colouring of G some pair at one of
two or three distances is monochromatic -- a weak property over a handful of
distances, which a multispindle can consume.  If it takes thirty, the graph is
nowhere near, and the count says by how much.

Ordered smallest class first, so each step adds as little as possible and the
count is as informative as it can be.  The last step, all 36 classes at once,
turns G's 7877 edges into 29221 -- average degree 37, where a graph is
6-chromatic for reasons that have nothing to do with the plane -- so the
interesting answer is a SMALL count, and a large one is nearly vacuous.  That
caveat is why the count is reported rather than the bare verdict.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def run(name, P, k):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    n = len(P)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    Eset = set(E)
    byd = defaultdict(list)
    for i in range(n - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in Eset:
                v = int(sq[off, 0])
                if v:
                    byd[Fr(v, D2)].append((i, j))
    clo = sorted((d for d in byd if closable_distance(d)),
                 key=lambda d: len(byd[d]))
    print(f"\n{name} at {k} colours: {n} pts, {len(E)} edges, {len(clo)} "
          f"closable classes, sizes {[len(byd[d]) for d in clo]}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    cum = list(E)
    for t, d in enumerate(clo, 1):
        cum = cum + byd[d]
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in cum:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        deg = 2 * len(cum) / n
        print(f"   {t:2d} classes (through D={d}): {len(cum)} edges, mean "
              f"degree {deg:.1f} -> "
              f"{'colours' if ok else '*** DOES NOT COLOUR ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if not ok:
            print(f"   >>> {t} classes suffice; the weak property holds over "
                  f"those distances", flush=True)
            return t, deg
    return None, None


run("G", build_G(K, as_graph=False), 5)
run("Sa", build_Sa(K), 4)
run("Sa", build_Sa(K), 5)
print("\nDONE", flush=True)
