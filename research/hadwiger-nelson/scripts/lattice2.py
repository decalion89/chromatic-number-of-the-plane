"""Out of de Grey's family: two-distance graphs on lattices.

Everything tested so far descends from one seed.  S gives Sa gives Y gives G,
and the closures, bites and thickenings of G are still G.  The family carries
the weak property at four colours and nothing at five, and the bite cannot
create one -- so the search has to leave the family, not grow it.

The mechanism does not require it.  The weak property is exactly "the graph on
this point set with BOTH distance 1 and distance d joined is not
k-colourable", and the point set is free.  If a lattice patch has a
6-chromatic two-distance graph while its unit-distance graph alone is
5-colourable, that lattice carries the weak property at five -- and lattices
are not de Grey's.

Single-distance lattices were ruled out here earlier and correctly: they are
always bipartite, since dx^2 + dy^2 = r forces dx + dy = r mod 2, so an odd r
flips the parity of i+j on every edge and an even one reduces to a scaled
copy.  That argument says nothing about TWO distances, where the two parity
classes interleave and the bipartition dies.

Triangular lattice first, since its squared distances are the Loeschian
numbers a^2 + ab + b^2 -- 1, 3, 4, 7, 9, 12, 13 -- and 3, 4, 7, 9 are all
closable over de Grey's own field: 4d^2-1 gives 11, 15, 27, 35.  Square
lattice after, where 2 and 4 are closable via 7 and 15.
"""
import sys, time, math
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()


def triangular(R):
    """Lattice points a*(1,0) + b*(1/2, sqrt3/2) inside radius R.  Squared
    distances are integers in the (a, b) form a^2 + ab + b^2."""
    pts = []
    n = int(R) + 2
    for a in range(-n, n + 1):
        for b in range(-n, n + 1):
            if a * a + a * b + b * b <= R * R:
                pts.append((a, b))
    return pts, lambda p, q: ((p[0] - q[0]) ** 2
                              + (p[0] - q[0]) * (p[1] - q[1])
                              + (p[1] - q[1]) ** 2)


def square(R):
    pts = []
    n = int(R) + 2
    for a in range(-n, n + 1):
        for b in range(-n, n + 1):
            if a * a + b * b <= R * R:
                pts.append((a, b))
    return pts, lambda p, q: (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2


def test(name, pts, d2f, D):
    n = len(pts)
    E, X = [], []
    for i in range(n - 1):
        for j in range(i + 1, n):
            v = d2f(pts[i], pts[j])
            if v == 1:
                E.append((i, j))
            elif v == D:
                X.append((i, j))
    if not X:
        return None
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E + X:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    verdict = ("colours" if ok else
               "*** DOES NOT COLOUR -- THE WEAK PROPERTY AT FIVE ***")
    print(f"  {name} R={len(pts)}pts  d^2={D}: {len(E)} unit + {len(X)} "
          f"second-distance edges -> {verdict}  [{time.time()-t0:.0f}s]",
          flush=True)
    return ok


for lat, gen in (("triangular", triangular), ("square", square)):
    print(f"\n=== {lat} lattice ===", flush=True)
    for R in (6, 9, 12):
        pts, d2f = gen(R)
        cand = sorted({d2f(pts[0], q) for q in pts} - {0, 1})
        clo = [d for d in cand if d <= 40 and closable_distance(Fr(d))]
        print(f" radius {R}: {len(pts)} points; closable second distances "
              f"{clo}  [{time.time()-t0:.0f}s]", flush=True)
        for D in clo:
            r = test(lat, pts, d2f, D)
            if r is False:
                print("  *** FOUND ***", flush=True)
                sys.exit(0)
print("\nDONE", flush=True)
