"""What a spindle actually glues.

The folk story is that spindling joins two copies "at the hinge": the points on
the radius-r circle, each to its own image.  On its own that story cannot be
right.  The hinge is a perfect matching between the two copies on R and the
identity elsewhere, so if the hinge were all there is, the glued graph would be
a subgraph of the Cartesian product H box K2 -- and by Sabidussi's theorem
chi(H box K2) = max(chi(H), 2) = chi(H).  A pure hinge can never raise the
chromatic number, no matter how many points are on the circle.

So the gain has to come from somewhere else: the points the two copies SHARE,
and the unit distances between them that are not hinge edges.  This measures
all three for de Grey's own two spindles.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Sb, build_Y, build_G, _rot_half_pi_pm
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph

F = Field((3, 5, 7, 11))
ORI = Point(F.zero(), F.zero())
PIV = Point(F.rational(-2), F.zero())

def report(name, A, B, centre, d2):
    sa, sb = set(A), set(B)
    shared = sa & sb
    hinge = sum(1 for p in A if (p - centre).norm2() == d2)
    onlyA = [p for p in A if p not in shared]
    onlyB = [p for p in B if p not in shared]
    cross = 0
    hinge_edges = 0
    rot = rotation_joining(d2, F).about(centre)
    for p in onlyA:
        img = rot(p)
        for q in onlyB:
            if p.is_unit_apart(q):
                cross += 1
                if q == img:
                    hinge_edges += 1
    print(f"{name}")
    print(f"   |A| = {len(A)}   |B| = {len(B)}   shared = {len(shared)}")
    print(f"   points on the radius-r circle of A: {hinge}")
    print(f"   cross edges between the private parts: {cross}")
    print(f"     of which hinge (p -- rho(p)): {hinge_edges}")
    print(f"     of which incidental: {cross - hinge_edges}"
          f"   ({100.0*(cross-hinge_edges)/max(cross,1):.1f}% of the interface)", flush=True)

Sa = build_Sa(F); Sb = build_Sb(F)
report("spindle 1:  Sa  u  rho_{r^2=4}(Sa)   (about the origin)", Sa, Sb, ORI, Fr(4))

Y = build_Y(F)
Ya = [_rot_half_pi_pm(F, +1).about(PIV)(p) for p in Y]
Yb = [_rot_half_pi_pm(F, -1).about(PIV)(p) for p in Y]
report("\nspindle 2:  Ya  u  Yb              (about (-2,0), r^2=16)", Ya, Yb, PIV, Fr(16))
