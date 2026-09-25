"""The menu of spindles available on a hub.

A rotation by 2*arcsin(1/(2r)) about a centre c sends every point at distance
exactly r from c to a point at distance exactly 1 from itself.  So gluing H to
its image under that rotation creates a *guaranteed* set of new unit edges --
one per point of H on the radius-r circle about c.  That is the hinge, and its
size is the only part of the gain that can be predicted before building.

de Grey's G is three such spindles chained:
    r^2 = 3   inside S itself    (needs sqrt 11)
    r^2 = 4   about the origin   (needs sqrt 15 = sqrt3 sqrt5)
    r^2 = 16  about (-2, 0)      (needs sqrt 63 = 3 sqrt7)
and Q(sqrt3,sqrt5,sqrt7,sqrt11) is exactly the field those three radii force.

This script lists, for a given hub, every (centre, radius) pair whose hinge is
non-trivial, together with the radical the spindle costs.  de Grey's own two
choices must appear in the list; the question is whether they are the best
ones in it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, pickle
from fractions import Fraction as Fr
from collections import Counter, defaultdict

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_S
from hn.geometry import Point, required_radical
from hn.field import Field

F = Field((3, 5, 7, 11))
which = sys.argv[1] if len(sys.argv) > 1 else "Sa"
pts = {"S": build_S, "Sa": build_Sa, "Y": build_Y}[which](F)
n = len(pts)
print(f"{which}: {n} points", flush=True)

def rat(e):
    """The rational value of a field element, or None if it is irrational."""
    if not e.is_rational():
        return None
    return Fr(e.c[0])

# candidate centres: the origin, every hub point, and de Grey's pivot (-2,0)
cents = {("origin", Point(F.zero(), F.zero()))}
cents.add(("pivot(-2,0)", Point(F.rational(-2), F.zero())))
for i, p in enumerate(pts):
    cents.add((f"v{i}", p))
cents = sorted(cents, key=lambda t: t[0])
print(f"{len(cents)} candidate centres", flush=True)

best = []
for name, c in cents:
    buck = Counter()
    for p in pts:
        d2 = rat((p - c).norm2())
        if d2 is None or d2 == 0:
            continue
        buck[d2] += 1
    for d2, m in buck.items():
        if m < 3:
            continue
        rad = required_radical(d2)
        best.append((m, d2, rad, name))

best.sort(key=lambda t: (-t[0], t[1]))
print(f"\n{len(best)} (centre,radius) pairs with hinge >= 3\n")
print(f"{'hinge':>6} {'r^2':>12} {'radical':>8}  {'in Q(3,11)?':>11}  centre")
seen = set()
shown = 0
for m, d2, rad, name in best:
    key = (m, d2, rad)
    if key in seen and shown > 60:
        continue
    seen.add(key)
    infield = "yes" if rad in (1, 3, 11, 33) else "NO (sqrt%d)" % rad
    star = "   <-- de Grey" if (d2 in (Fr(3), Fr(4), Fr(16)) and name in ("origin", "pivot(-2,0)")) else ""
    print(f"{m:>6} {str(d2):>12} {rad:>8}  {infield:>11}  {name}{star}")
    shown += 1
    if shown >= 70:
        break

# the aggregate: which radii are rich anywhere at all
agg = defaultdict(int)
for m, d2, rad, name in best:
    agg[(d2, rad)] = max(agg[(d2, rad)], m)
print(f"\n--- best hinge per radius, radicals already in Q(sqrt3,sqrt11) ---")
free = sorted(((m, d2, rad) for (d2, rad), m in agg.items() if rad in (1, 3, 11, 33)),
              reverse=True)[:20]
for m, d2, rad in free:
    print(f"  hinge {m:>4}  r^2 = {str(d2):<12} radical {rad}")
print(f"\n--- best hinge per radius, radicals OUTSIDE it (what they would cost) ---")
paid = sorted(((m, d2, rad) for (d2, rad), m in agg.items() if rad not in (1, 3, 11, 33)),
              reverse=True)[:20]
for m, d2, rad in paid:
    print(f"  hinge {m:>4}  r^2 = {str(d2):<12} radical {rad}")
