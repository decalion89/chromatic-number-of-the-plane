"""The template at five colours, with the centre where the template puts it.

Verified above, and every piece of it shares ONE centre, the origin:

  Sa  = the D6 closure of S about the origin        (397 pts, 4-chromatic)
  rho = the bite for the ring D=4, about the origin (cos 7/8, sin sqrt15/8)
  Y   = Sa u rho(Sa), less two points               (791 pts)
        -> forces the ring's ANTIPODAL pair (-2,0),(2,0), squared distance 16
  G   = Y u sigma(Y), sigma the D=16 spindle about ONE END, (-2,0)
        (1581 pts, 5-chromatic)

Measured, not assumed: Y forces that pair at four colours (UNSAT in 210s),
Sa alone does not, and Y does not force it at five.

G* as built earlier is closed about (-2,0), which is the spindle pivot, not the
closure centre -- the wrong point for this.  So close about a centre and use
THAT centre for the ring, the bite and the antipodal map, all three.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Y
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.geometry import _rot60
from hn.graph import build_graph
from hn.homcol import doubly_usable_ring
from pysat.solvers import Solver

k = 5
t0 = time.time()
ORIGIN = Point(F.zero(), F.zero())
PIV = Point(F.rational(-2), F.zero())


def dihedral(pts, about):
    rot = _rot60(F).about(about)
    out, seen = [], set()
    for refl in (False, True):
        for j in range(6):
            for p in pts:
                q = p
                for _ in range(j):
                    q = rot(q)
                if refl:
                    q = Point(q.x, about.y + (about.y - q.y))
                if q not in seen:
                    seen.add(q)
                    out.append(q)
    return out


def rings(P, C):
    out = defaultdict(list)
    for i, p in enumerate(P):
        d = p.dist2(C)
        if all(x == 0 for x in d.c[1:]) and d.c[0] != 0:
            out[d.c[0]].append(i)
    return out


def study(name, B, C):
    R = rings(B, C)
    du = sorted(d for d in R if doubly_usable_ring(d))
    print(f"\n=== {name}: {len(B)} pts, {len(R)} rational rings about the "
          f"centre, {len(du)} doubly usable: {du}  [{time.time()-t0:.0f}s]",
          flush=True)
    for D in du:
        rot = rotation_joining(D, F).about(C)
        seen, U = set(), []
        for p in list(B) + [rot(p) for p in B]:
            if p not in seen:
                seen.add(p)
                U.append(p)
        g = build_graph(U)
        E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
        cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        print(f"  {name} D={D}: union {g.n} pts, {len(E)} edges, solving"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        sv = Solver(name="cd15", bootstrap_with=cls)
        if not sv.solve():
            print(f"  *** {name} D={D}: UNION IS NOT 5-COLOURABLE ***",
                  flush=True)
            return
        idx = {p: i for i, p in enumerate(U)}
        anti = []
        for i in R[D]:
            p = B[i]
            q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
            if q in idx and idx[p] < idx[q]:
                anti.append((idx[p], idx[q]))
        print(f"    {len(anti)} antipodal pairs on the ring (squared "
              f"distance {4*D}); testing  [{time.time()-t0:.0f}s]", flush=True)
        hits = []
        for i, j in anti:
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)]):
                hits.append((i, j))
                print(f"  *** FORCED SAME: {i},{j} at D'={4*D} ***",
                      flush=True)
        print(f"    {len(hits)} of {len(anti)} forced"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        if hits:
            return


G = build_G(F, as_graph=False)
study("G closed about the origin", dihedral(G, ORIGIN), ORIGIN)
study("G closed about (-2,0)", dihedral(G, PIV), PIV)
print("\nDONE", flush=True)
