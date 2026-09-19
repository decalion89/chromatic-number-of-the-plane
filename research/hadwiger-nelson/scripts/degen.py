"""Why chi(confined) will not move: measure the degeneracy that caps it.

Adding hexagons takes the auxiliaries from 126 to 798 and leaves chi(confined)
at [2,3] throughout, which asks for an explanation rather than more search.  A
graph of degeneracy d is (d+1)-colourable, so if the confined set is
2-degenerate its chromatic number cannot exceed 3 no matter how large it gets
-- and 4 is what pressure 3 at five colours needs.

Degeneracy is the right statistic here for a second reason: Erdos-Rubin-Taylor
says a 2-degenerate graph is also 3-CHOOSABLE, so the lists of size 3 would
complete even if the colours were not all the same.  Both halves of the
obstruction would then be the same number.
"""
import sys, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn.field import Field
from hn.geometry import Point, Rotation

FLD = Field((3, 11))
ONE, ZERO = FLD.rational(1), FLD.zero()
ROT60 = Rotation(FLD.rational(Fraction(1, 2)),
                 FLD.sqrt(3) * FLD.rational(Fraction(1, 2)))
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


offsets, q = [], Point(ONE, ZERO)
for _ in range(8):
    offsets.append(q)
    q = rot(CH, q)

print(f"{'hexagons':>9} {'aux':>5} {'confined':>9} {'max deg':>8} "
      f"{'degeneracy':>11}  cap on chi")
print("-" * 58)
for t in range(3, 8):
    hexes = [hexagon(s) for s in offsets[:t]]
    aux = {}
    for a in range(t):
        for b in range(a, t):
            for i, u in enumerate(hexes[a]):
                for j, v in enumerate(hexes[b]):
                    if (a, i) >= (b, j):
                        continue
                    s = Point(u.x + v.x, u.y + v.y)
                    if s.norm2() in (ONE, ZERO):
                        continue
                    aux.setdefault((s.x, s.y), set()).add((a, i))
                    aux[(s.x, s.y)].add((b, j))
    pts = list(aux)
    adj = {i: set() for i in range(len(pts))}
    for i in range(len(pts)):
        pi = Point(*pts[i])
        for j in range(i + 1, len(pts)):
            if pi.is_unit_apart(Point(*pts[j])):
                adj[i].add(j)
                adj[j].add(i)

    worstdeg, worstsize, worstmax = 0, 0, 0
    for bits in itertools.product((0, 1), repeat=t):
        col = {(a, i): (i + bits[a]) % 2 for a in range(t) for i in range(6)}
        conf = {i for i in range(len(pts))
                if {col[n] for n in aux[pts[i]]} == {0, 1}}
        d = degeneracy(adj, conf)
        if d > worstdeg:
            worstdeg, worstsize = d, len(conf)
            worstmax = max((len(adj[v] & conf) for v in conf), default=0)
    print(f"{t:>9} {len(pts):>5} {worstsize:>9} {worstmax:>8} "
          f"{worstdeg:>11}  chi <= {worstdeg + 1}", flush=True)
