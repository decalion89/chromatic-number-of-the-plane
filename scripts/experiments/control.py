"""Positive control for the forcing test, on the smallest thing that forces.

A rhombus of two unit triangles -- A, two shoulders, and the far tip D at
distance sqrt(3) -- forces A and D to share a colour at three colours: the
triangle A,B,C uses all three, the triangle B,C,D then leaves D only the one
A already has.  That is the whole engine of the Moser spindle, and if the
census code cannot see it the zero it reports elsewhere is worthless.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as K, Point
from hn.graph import build_graph
from hn.homcol import forced_same, forced_different
from pysat.solvers import Solver

half, r3, r11 = K.rational(Fr(1, 2)), K.sqrt(3), K.sqrt(11)
sixth = K.rational(Fr(1, 6))


def rot(p, c, s):
    return Point(p.x * c - p.y * s, p.x * s + p.y * c)


A = Point(K.zero(), K.zero())
D = Point(r3, K.zero())
B = Point(r3 * half, half)
C = Point(r3 * half, -half)
rh = [A, B, C, D]
# second rhombus, hinged at A, turned until the tips are a unit apart:
# |D|=|D'|=sqrt3 and |D-D'|=1 needs cos = 5/6, sin = sqrt11/6.
c, s = K.rational(Fr(5, 6)), r11 * sixth
sp = rh + [rot(p, c, s) for p in (B, C, D)]


def report(name, pts, u, v):
    g = build_graph(pts)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    print(f"{name}: {g.n} points, {len(E)} edges, "
          f"d^2({u},{v}) = {pts[u].dist2(pts[v])}")
    for k in (3, 4):
        cls = [[1 + w * k + cc for cc in range(k)] for w in range(g.n)]
        for a, b in E:
            for cc in range(k):
                cls.append([-(1 + a * k + cc), -(1 + b * k + cc)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        fs = forced_same(sv, u, v, k) if ok else None
        fd = forced_different(sv, u, v, k) if ok else None
        sv.delete()
        print(f"   k={k}: colourable={ok}  forced same={fs}  forced diff={fd}")


report("rhombus", rh, 0, 3)
report("spindle", sp, 3, 6)
