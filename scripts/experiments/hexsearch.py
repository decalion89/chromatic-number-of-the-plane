"""How far the confined set gets, as a function of how many hexagons there are.

The quantity that gates pressure 3 at five colours is now known exactly: the
CONFINED set -- the auxiliaries whose two circle neighbours are coloured
differently -- must be 4-chromatic in EVERY proper 2-colouring of the circle.
de Grey's Sa reaches a maximum of 3 and a minimum of 2 over its 32
orientations, so the gap is real and this is a graded objective rather than
the binary refuses/does-not that every earlier search used.

The circle is a union of hexagons, one per angle offset, each a 6-cycle with
exactly two proper 2-colourings; the auxiliaries are all Minkowski sums u + v
of circle points in different hexagons, since a point one away from u and v is
the pivot reflected across the chord uv.  So a configuration is just a list of
angle offsets, and the objective is

    min over orientations of chi(confined set)

which has to reach 4.  This measures it as hexagons are added, over the angle
family that de Grey's construction generates and over random ones, to see
whether the number moves at all.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn.field import Field
from hn.geometry import Point, Rotation
from pysat.solvers import Solver

FLD = Field((3, 11))
ONE, ZERO = FLD.rational(1), FLD.zero()
ROT60 = Rotation(FLD.rational(Fraction(1, 2)),
                 FLD.sqrt(3) * FLD.rational(Fraction(1, 2)))

# cos/sin of the offsets de Grey's construction generates, in exact arithmetic
CT = (FLD.rational(Fraction(5, 6)), FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
CH = (FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
      FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))


def rot(cs, p):
    c, s = cs
    return Point(p.x * c - p.y * s, p.x * s + p.y * c)


def hexagon(seed):
    q, out = seed, []
    for _ in range(6):
        out.append(q)
        q = ROT60(q)
    return out


def chi(adj, verts, hi=5):
    vs = sorted(verts)
    if not vs:
        return 0
    pos = {v: i for i, v in enumerate(vs)}
    es = [(pos[a], pos[b]) for a in vs for b in adj.get(a, ()) if b in pos and a < b]
    for k in range(1, hi + 1):
        cls = [[1 + i * k + c for c in range(k)] for i in range(len(vs))]
        for a, b in es:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return hi + 1


def evaluate(seeds, label):
    hexes = [hexagon(s) for s in seeds]
    circle = [p for h in hexes for p in h]
    if len({(p.x, p.y) for p in circle}) < len(circle):
        return None                      # coincident hexagons are not a design
    # A point one away from circle points u and v is the pivot reflected
    # across the chord uv, which is u + v: |u + v - u| = |v| = 1 for ANY pair
    # on the circle, not just for chords of length one.  So Minkowski summing
    # the circle with itself enumerates every auxiliary there is, and pairs
    # WITHIN one hexagon count -- an earlier version took only cross-hexagon
    # pairs and dropped the always-confined points at sqrt 3 from the pivot.
    aux = {}
    for a in range(len(hexes)):
        for b in range(a, len(hexes)):
            for i, u in enumerate(hexes[a]):
                for j, v in enumerate(hexes[b]):
                    if (a, i) >= (b, j):
                        continue
                    q = Point(u.x + v.x, u.y + v.y)
                    if q.norm2() == ONE or q.norm2() == ZERO:
                        continue   # on the circle, or the pivot itself (v = -u)
                    aux.setdefault((q.x, q.y), set()).add((a, i))
                    aux[(q.x, q.y)].add((b, j))
    pts = list(aux)
    adj = {i: set() for i in range(len(pts))}
    for i in range(len(pts)):
        pi = Point(*pts[i])
        for j in range(i + 1, len(pts)):
            if pi.is_unit_apart(Point(*pts[j])):
                adj[i].add(j)
                adj[j].add(i)

    worst, best = 99, 0
    for bits in itertools.product((0, 1), repeat=len(hexes)):
        col = {(a, i): (i + bits[a]) % 2 for a in range(len(hexes))
               for i in range(6)}
        confined = {i for i in range(len(pts))
                    if {col[n] for n in aux[pts[i]]} == {0, 1}}
        c = chi(adj, confined)
        worst, best = min(worst, c), max(best, c)
    print(f"  {label}: {len(hexes)} hexagons, {len(pts)} auxiliaries, "
          f"chi(confined) in [{worst}, {best}]"
          + ("   *** MIN 4 ***" if worst >= 4 else ""), flush=True)
    return worst


print("de Grey's own angle family, as hexagons are added:")
base = [Point(ONE, ZERO), Point(*CH), Point(*CT)]
for t in range(2, 4):
    evaluate(base[:t], f"offsets[0:{t}]")

print("\nrandom extra offsets on top of the three (exact rotations of them):")
rng = random.Random(5)
extra = [rot(CH, rot(CH, Point(ONE, ZERO))),        # theta
         rot(CT, Point(*CH)),                        # 3 theta / 2
         rot(CT, rot(CT, Point(ONE, ZERO))),         # 2 theta
         rot(CH, rot(CT, Point(*CT)))]               # 5 theta / 2
for k in range(1, len(extra) + 1):
    evaluate(base + extra[:k], f"three + {k}")
