"""sigma: the exact map that lands a point on its own 1/sqrt3 ring, a unit away.

The cross angle between the radius-1 ring and the radius-1/sqrt3 ring is
theta0 = arccos(sqrt3/6), with sin theta0 = sqrt(11/12) = sqrt33/6.  Both live
in Q(sqrt3, sqrt11), so the rotation by theta0 composed with the scaling by
1/sqrt3 is an exact map of the field this project already works in:

    sigma(x, y) = ( (x - sqrt11 y)/6 , (sqrt11 x + y)/6 )

and two lines of algebra give everything.  |sigma(u)|^2 = (x^2+y^2)/3, so sigma
lands every point on the circle of 1/sqrt3 times its radius; and
|u - sigma(u)|^2 = x^2 + y^2, so sigma moves each point by exactly its own
radius -- a unit, for a point of the unit circle.  The conjugate map sigmabar,
with -sqrt11, does the same on the other side.

So the cross edges that the triangular lattice cannot supply -- its angles are
multiples of 30 degrees and theta0 is not one -- are constructible after all,
and this is why the centre-closure failed: it manufactured 1/sqrt3 pairs by the
thousand but never a single edge between the two rings of the same hub.

The smallest test: h at the origin, the unit hexagon N around it, and
sigma(N) with sigmabar(N) on the 1/sqrt3 ring.  Nineteen points.  Each u in N
then has exactly two ring neighbours, and if those two differ in colour u is
left a single choice out of the ring's three -- which is the first time in this
project that a neighbourhood point has been squeezed at all.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from itertools import combinations
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
F = Field((3, 11))
half = F.rational(Fr(1, 2)); r3 = F.sqrt(3); r11 = F.sqrt(11)
sixth = F.rational(Fr(1, 6))
def sigma(p, conj=False):
    s = -r11 if conj else r11
    return Point((p.x - s * p.y) * sixth, (s * p.x + p.y) * sixth)
O = Point(F.rational(0), F.rational(0))
hexa = []
cx, cy = half, r3 * half
x, y = F.rational(1), F.rational(0)
for _ in range(6):
    hexa.append(Point(x, y))
    x, y = x * cx - y * cy, x * cy + y * cx
one = F.rational(1); third = F.rational(Fr(1, 3))
assert all((p - O).norm2() == one for p in hexa)
assert all((p - sigma(p)).norm2() == one for p in hexa), "sigma must move by one"
assert all((sigma(p) - O).norm2() == third for p in hexa), "sigma must land on 1/sqrt3"
print(f"sigma verified exactly on the hexagon   [{time.time()-t0:.0f}s]", flush=True)

pts = {}
for p in [O] + hexa + [sigma(p) for p in hexa] + [sigma(p, True) for p in hexa]:
    pts[(round(p.fx, 9), round(p.fy, 9))] = p
G = build_graph(list(pts.values())); n = G.n
E = list(G.edges())
h = next(i for i in range(n) if (G.vertices[i] - O).norm2() == F.rational(0))
ring = [i for i in range(n) if (G.vertices[i] - O).norm2() == third]
nb = [i for i in range(n) if (G.vertices[i] - O).norm2() == one]
print(f"  n={n} edges={len(E)}; |N(h)|={len(nb)} |ring|={len(ring)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def chi(verts, edges, cap=6):
    for k in range(1, cap + 1):
        idx = {v: i for i, v in enumerate(verts)}
        X = lambda v, c: 1 + idx[v] * k + c
        cnf = [[X(v, c) for c in range(k)] for v in verts]
        for a, b in edges:
            for c in range(k):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(); s.delete()
        if ok: return k
    return cap + 1

S = set(nb) | set(ring)
Eloc = [(a, b) for a, b in E if a in S and b in S]
print(f"  chi(N(h) u ring) = {chi(sorted(S), Eloc)}  on {len(S)} points, "
      f"{len(Eloc)} edges   [{time.time()-t0:.0f}s]", flush=True)
print(f"  chi(whole configuration) = {chi(list(range(n)), E)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

K = 5
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for a, b in E:
    for c in range(K):
        cnf.append([-X(a, c), -X(b, c)])
for v in ring:
    for c in range(K):
        cnf.append([-X(h, c), -X(v, c)])
s = Solver(name="cd19", bootstrap_with=cnf)
print(f"  h can avoid its whole ring at five colours: {s.solve()}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
s.delete()
