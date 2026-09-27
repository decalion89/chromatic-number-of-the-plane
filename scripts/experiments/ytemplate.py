"""Ask de Grey's forcer which pair it forces.  It has never been asked.

Everything here has guessed at what a forcer looks like -- ring population,
which radicals a distance pays, how often a pair agrees in sampled colourings.
Y is a forcer, it is 791 points, and the exact test now costs a millisecond,
so the guessing can stop: enumerate every non-edge pair of Y at a rational
squared distance and test each one exactly at FOUR colours.

The pair (i, j) is forced-same iff adding the edge (i, j) makes the graph
uncolourable, which by colour symmetry is one call with assumptions
[x_{i,0}, -x_{j,0}] coming back UNSAT.

Whatever comes out is the template: its distance, its position relative to the
pivot (-2, 0) that build_G rotates about, and how many such pairs there are.
The same scan at FIVE colours on the same graph is the control -- Y is
5-colourable with room to spare, so it should force nothing there.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def scan(name, P, k, only_closable=True):
    g = build_graph(P)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    pairs = []
    for i in range(n - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in E:
                pairs.append((i, j, Fr(int(sq[off, 0]), D2)))
    keep = [p for p in pairs if not only_closable or closable_distance(p[2])]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"{name}: NOT {k}-colourable", flush=True)
        return n, []
    forced = [(i, j, d) for i, j, d in keep
              if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    print(f"{name} at {k} colours: {n} pts, {len(E)} edges, "
          f"{len(pairs)} rational non-edge pairs, {len(keep)} tested, "
          f"{len(forced)} FORCED  [{time.time()-t0:.0f}s]", flush=True)
    return n, forced


# --- the template: Y at four colours ---------------------------------------
Y = build_Y(F)
n, fY = scan("Y", Y, 4)
if fY:
    byd = Counter(d for _, _, d in fY)
    print(f"  distances: {dict(byd)}", flush=True)
    piv = Point(F.rational(-2), F.zero())
    for i, j, d in fY[:12]:
        a, b = Y[i], Y[j]
        print(f"    {i},{j} at D={d}: "
              f"({a.fx:+.4f},{a.fy:+.4f}) and ({b.fx:+.4f},{b.fy:+.4f}); "
              f"|a-piv|^2={a.dist2(piv)}, |b-piv|^2={b.dist2(piv)}", flush=True)

# --- the same scan without the closability filter ---------------------------
print("", flush=True)
scan("Y (all rational, not just closable)", Y, 4, only_closable=False)

# --- the control: Y at five colours, and Sa at four ------------------------
print("", flush=True)
scan("Y", Y, 5, only_closable=False)
print("", flush=True)
scan("Sa", build_Sa(F), 4, only_closable=False)
print("\nDONE", flush=True)
