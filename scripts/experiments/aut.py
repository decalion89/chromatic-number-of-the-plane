"""The isometries a point set actually admits.

Symmetrising the 951-vertex graph under D6 -- Sa's group -- loosened it, which
is what you expect if D6 is not ITS group.  So ask what group it does carry.

Any isometry of a finite planar set fixes the centroid, so every candidate is a
rotation or a reflection about that point, and each is determined by where it
sends one chosen vertex.  Division is avoided throughout: the map is applied in
the SCALED form r -> (co*x -+ si*y, si*x +- co*y) with co, si taken before
dividing by |a|^2, and compared against the set of |a|^2 * (q - centroid).
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point
from hn.field import Field

SRC = sys.argv[1]
if SRC == "Sa":
    F = Field((3, 11, 247)); P = build_Sa(F)
else:
    d = json.load(open(f"{HN_DIR}/data/{SRC}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
n = len(P)
sx, sy = P[0].x, P[0].y
for p in P[1:]:
    sx, sy = sx + p.x, sy + p.y
inv = F.rational(Fr(1, n))
C = Point(sx * inv, sy * inv)
rel = [p - C for p in P]
print(f"{SRC}: {n} points, centroid {C.approx()}", flush=True)

anchor = max(rel, key=lambda r: float(r.x * r.x + r.y * r.y))
na = anchor.x * anchor.x + anchor.y * anchor.y
scaled = {Point(r.x * na, r.y * na) for r in rel}
t0 = time.time()
rot, ref = 0, 0
for b in rel:
    if b.x * b.x + b.y * b.y != na:
        continue
    ax, ay, bx, by = anchor.x, anchor.y, b.x, b.y
    for sign in (1, -1):
        if sign > 0:
            co = ax * bx + ay * by; si = ax * by - ay * bx
        else:
            co = ax * bx - ay * by; si = ax * by + ay * bx
        ok = True
        for r in rel:
            if sign > 0:
                q = Point(co * r.x - si * r.y, si * r.x + co * r.y)
            else:
                q = Point(co * r.x + si * r.y, si * r.x - co * r.y)
            if q not in scaled:
                ok = False; break
        if ok:
            if sign > 0: rot += 1
            else: ref += 1
print(f"  isometry group: {rot} rotations (including the identity), "
      f"{ref} reflections -> order {rot + ref}   [{time.time()-t0:.0f}s]", flush=True)
