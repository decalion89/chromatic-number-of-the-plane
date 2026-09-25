"""de Grey's template, run at five colours -- exactly, no sampling.

The template, now established rather than guessed:

  Sa is the 12-fold dihedral closure of S about the origin.  It carries a ring
  of squared radius D = 4.  Rotating Sa about the ORIGIN by the angle that
  makes that ring bite itself -- cos 7/8, sin sqrt(15)/8, which needs
  sqrt(4D - 1) = sqrt(15) -- gives Sb, and Y = Sa u Sb forces the ring's
  ANTIPODAL pair (-2, 0), (2, 0) monochromatic in every 4-colouring.  That
  pair is at squared distance 4D = 16, and spindling it needs
  sqrt(16D - 1) = sqrt(63) = 3 sqrt(7).  Both radicals are in the field, which
  is why D = 4 and only D = 4 works: the ring pays twice.

So the object to look for at five colours is not "a forced pair somewhere".
It is the antipodal pair of a doubly-usable ring in B u rho_D(B), and testing
it is ONE SAT call: adding the edge makes the union uncolourable iff the pair
is forced.  That is cheap enough to sweep every base and every ring it has.

Each union also gets the full exhaustive treatment while it is in memory --
every non-edge pair at a rational closable distance, tested exactly -- so a
forced pair that is not the antipodal one is not missed either.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.fast import IntBasis
from hn.homcol import closable_distance, doubly_usable_ring
from pysat.solvers import Solver

k = 5
t0 = time.time()
O = Point(F.zero(), F.zero())


def rings(P):
    """Rational squared radii about the origin, with their point lists."""
    out = defaultdict(list)
    for i, p in enumerate(P):
        d = p.norm2()
        if all(x == 0 for x in d.c[1:]) and d.c[0] != 0:
            out[d.c[0]].append(i)
    return out


def closable_pairs(P, E):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    dim, D2 = basis.dim, basis.D * basis.D
    out = []
    for i in range(len(P) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) in E:
                continue
            dd = Fr(int(sq[off, 0]), D2)
            if closable_distance(dd):
                out.append((i, j, dd))
    return out


def solver_for(P):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    return Solver(name="cd15", bootstrap_with=cls), E, g.n


def study(name, B):
    R = rings(B)
    du = sorted(d for d in R if doubly_usable_ring(d))
    print(f"\n=== {name}: {len(B)} pts, {len(R)} rational rings, "
          f"{len(du)} doubly usable: {du}  [{time.time()-t0:.0f}s]",
          flush=True)
    for D in du:
        rot = rotation_joining(D, F)          # the biting rotation about O
        seen, U = set(), []
        for p in list(B) + [rot(p) for p in B]:
            if p not in seen:
                seen.add(p)
                U.append(p)
        sv, E, n = solver_for(U)
        if not sv.solve():
            print(f"  *** {name} D={D}: UNION IS NOT 5-COLOURABLE "
                  f"({n} pts) ***", flush=True)
            return
        idx = {p: i for i, p in enumerate(U)}
        # the ring's antipodal pairs, in the union
        anti = []
        for i in R[D]:
            p = B[i]
            q = Point(-p.x, -p.y)
            if q in idx and idx[p] < idx[q]:
                anti.append((idx[p], idx[q]))
        hits = [(i, j) for i, j in anti
                if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
        print(f"  {name} D={D}: union {n} pts, {len(E)} edges, "
              f"{len(anti)} antipodal pairs on the ring, {len(hits)} FORCED"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if hits:
            print(f"  *** FORCED ANTIPODAL PAIR: {hits[:5]} at D'={4*D} ***",
                  flush=True)
            return
        if n <= 6000:
            cp = closable_pairs(U, E)
            f2 = [(i, j, d) for i, j, d in cp
                  if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
            print(f"    exhaustive: {len(cp)} closable non-edge pairs, "
                  f"{len(f2)} FORCED  [{time.time()-t0:.0f}s]", flush=True)
            if f2:
                print(f"  *** FORCED: {f2[:5]} ***", flush=True)
                return
        sv.delete()


study("Y", build_Y(F))
study("G", build_G(F, as_graph=False))
print("\nDONE", flush=True)
