"""The unit square: a width-two hub disjunction consumed by TWO copies.

The pigeonhole closure needs m > w copies and then, because the adversary picks
which pair repeats, 2r sin(k psi/2) = 1 for every k at once -- which forces
m = 3, psi = 120 degrees, r = 1/sqrt3.  But pigeonhole is not the only way to
close a disjunction.  With m = 2 and w = 2 one can instead demand that ALL FOUR
cross pairs clash, and that has a solution.

Write v1, v2 for the partners, rho for the rotation by phi about h.  Then
|v_s - rho v_s| = 1 forces r1 = r2 = r with 2r sin(phi/2) = 1, and the two cross
conditions force cos(a2 - a1 + phi) = cos(phi - a2 + a1), so a2 - a1 = 180
degrees: the partners are ANTIPODAL about h.  Substituting, 1 + cos phi =
1 - cos phi, so phi = 90 and r = 1/sqrt2.

Which is the unit square.  h is its centre, v1 and v2 are one diagonal, and
rho v1, rho v2 are the other; every point of one diagonal is a unit from every
point of the other, because those are the four sides.  So

    if some vertex h of a unit-distance graph G is the centre of a unit square
    whose diagonal {v1, v2} lies in G, and every 5-colouring of G gives
    c(h) = c(v1) or c(h) = c(v2), then G u rho_90(G) has no 5-colouring.

Two copies.  And the 90-degree rotation is rational, so the field never grows --
unlike the 120-degree version, which needs sqrt3.

What it needs instead is a pair at distance sqrt2, and d^2 = 2 is not a
Loeschian number: a^2 + ab + b^2 never equals 2, since 2 is inert in the
Eisenstein integers.  Every carrier in this project lives in a triangular
lattice, so none of them contains a single such pair -- measured, zero 90-degree
pairs in all four graphs scanned.  This construction needs a carrier the
project has never built.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point

t0 = time.time()
F = Field((3,))
def P(a, b): return Point(F.rational(Fr(a)), F.rational(Fr(b)))
h = Point(F.rational(Fr(1, 2)), F.rational(Fr(1, 2)))
v1, v2 = P(0, 0), P(1, 1)
def rot90(p):                      # about h: (x,y) -> (1 - y, x), exactly rational
    return Point(F.rational(1) - p.y, p.x)
r1, r2 = rot90(v1), rot90(v2)
one = F.rational(1); half = F.rational(Fr(1, 2)); two = F.rational(2)
print("the unit square, checked exactly:", flush=True)
print(f"  h  = ({h.fx},{h.fy})   v1 = ({v1.fx},{v1.fy})   v2 = ({v2.fx},{v2.fy})",
      flush=True)
print(f"  rho v1 = ({r1.fx},{r1.fy})   rho v2 = ({r2.fx},{r2.fy})", flush=True)
assert rot90(h) == h or (rot90(h) - h).norm2() == F.rational(0), "rho must fix h"
for nm, a, b in (("v1-rho v1", v1, r1), ("v1-rho v2", v1, r2),
                 ("v2-rho v1", v2, r1), ("v2-rho v2", v2, r2)):
    d2 = (a - b).norm2()
    print(f"  |{nm}|^2 = {d2.fx if hasattr(d2,'fx') else float(d2)}  "
          f"-> {'UNIT EDGE' if d2 == one else 'NOT A UNIT EDGE'}", flush=True)
    assert d2 == one
print(f"  |h - v1|^2 = {float((h-v1).norm2())} = 1/2, so r = 1/sqrt2", flush=True)
print(f"  |v1 - v2|^2 = {float((v1-v2).norm2())} = 2, so the diagonal is sqrt2",
      flush=True)
assert (h - v1).norm2() == half and (h - v2).norm2() == half
assert (v1 - v2).norm2() == two
print(f"\n  all four cross pairs are unit edges: the closure holds with TWO "
      f"copies   [{time.time()-t0:.0f}s]", flush=True)
print("  and 2 is not Loeschian (a^2+ab+b^2 = 2 has no integer solution), so no"
      "\n  triangular-lattice carrier contains the diagonal this needs.", flush=True)
for a in range(-6, 7):
    for b in range(-6, 7):
        assert a*a + a*b + b*b != 2
print("  checked over a^2+ab+b^2 for |a|,|b| <= 6.", flush=True)
