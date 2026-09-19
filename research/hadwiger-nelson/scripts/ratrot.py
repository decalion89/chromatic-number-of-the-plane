"""Is degeneracy 2 a property of de Grey's angles, or of hexagons as such?

The confined set is 2-degenerate throughout de Grey's own family, hexagons at
multiples of theta/2 where cos theta = 5/6.  That family is one curve through
a large space, so the question is whether any family does better -- degeneracy
3 is what the list argument needs, and it is the only thing that needs
finding.

Pythagorean rotations give a completely different family: cos and sin both
RATIONAL, so the whole configuration lives in Q(sqrt 3) -- hexagons need only
zeta_6 -- rather than in a quadratic extension.  The angles are unrelated to
the Moser rotation and to each other, which is exactly the point.
"""
import sys, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn.field import Field
from hn.geometry import Point, Rotation

FLD = Field((3,))
ONE, ZERO = FLD.rational(1), FLD.zero()
ROT60 = Rotation(FLD.rational(Fraction(1, 2)),
                 FLD.sqrt(3) * FLD.rational(Fraction(1, 2)))

PYTH = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29),
        (9, 40, 41), (12, 35, 37), (28, 45, 53)]


def rotation(a, b, c):
    return Rotation(FLD.rational(Fraction(a, c)), FLD.rational(Fraction(b, c)))


def hexagon(seed):
    q, out = seed, []
    for _ in range(6):
        out.append(q)
        q = ROT60(q)
    return out


def degeneracy(adj, verts):
    rem = {v: len(adj[v] & verts) for v in verts}
    d = 0
    while rem:
        v = min(rem, key=rem.get)
        d = max(d, rem[v])
        for u in adj[v]:
            if u in rem:
                rem[u] -= 1
        del rem[v]
    return d


def evaluate(seeds, label):
    hexes = [hexagon(s) for s in seeds]
    circle = [p for h in hexes for p in h]
    if len({(p.x, p.y) for p in circle}) < len(circle):
        print(f"  {label}: coincident hexagons, skipped", flush=True)
        return
    aux = {}
    for a in range(len(hexes)):
        for b in range(a, len(hexes)):
            for i, u in enumerate(hexes[a]):
                for j, v in enumerate(hexes[b]):
                    if (a, i) >= (b, j):
                        continue
                    q = Point(u.x + v.x, u.y + v.y)
                    if q.norm2() in (ONE, ZERO):
                        continue
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
    best = 0
    for bits in itertools.product((0, 1), repeat=len(hexes)):
        col = {(a, i): (i + bits[a]) % 2 for a in range(len(hexes))
               for i in range(6)}
        conf = {i for i in range(len(pts))
                if {col[n] for n in aux[pts[i]]} == {0, 1}}
        best = max(best, degeneracy(adj, conf))
    print(f"  {label}: {len(hexes)} hexagons, {len(pts)} auxiliaries, "
          f"best degeneracy {best}"
          + ("   *** 3 OR MORE ***" if best >= 3 else ""), flush=True)


print("powers of one Pythagorean rotation:")
for a, b, c in PYTH[:4]:
    r = rotation(a, b, c)
    seeds, q = [], Point(ONE, ZERO)
    for _ in range(5):
        seeds.append(q)
        q = r(q)
    evaluate(seeds, f"({a},{b},{c})^k, k<5")

print("\nmixed Pythagorean rotations, one hexagon each:")
for k in (3, 4, 5):
    for combo in itertools.combinations(PYTH, k - 1):
        seeds = [Point(ONE, ZERO)]
        for a, b, c in combo:
            seeds.append(rotation(a, b, c)(Point(ONE, ZERO)))
        evaluate(seeds, "+".join(f"{a}/{c}" for a, b, c in combo))
    if k == 4:
        break
