"""The confined set is indexed by a CUT, and that is a strong restriction.

Colour hexagon a with parity bits[a]: the circle point at position i in it
takes colour (i + bits[a]) mod 2.  An auxiliary is the sum u + v of circle
points at positions i in hexagon a and j in hexagon b, and it is confined
exactly when

    i + j + bits[a] + bits[b]  is odd.

Two consequences, and both bound what any orientation can do.

Within ONE hexagon the bits cancel, so confinement depends on i + j alone: the
same-hexagon auxiliaries are confined in EVERY orientation.  Those are the
points at sqrt 3 from the pivot, and their own graph is a union of paths --
bipartite, so they contribute at most 2 to any chromatic number.

Across hexagons only eps_ab = bits[a] XOR bits[b] matters, and eps is a CUT of
the complete graph K_t: eps_ab + eps_bc + eps_ac = 0 for every triple.  So the
2^t orientations give only 2^(t-1) distinct confined sets, and they are not
arbitrary -- no orientation can pick the pair-parities independently.

This checks both claims on the built configurations rather than asserting them.
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
CT = (FLD.rational(Fraction(5, 6)), FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
CH = (FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
      FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))


def hexagon(seed):
    q, out = seed, []
    for _ in range(6):
        out.append(q)
        q = ROT60(q)
    return out


t = 3
hexes = [hexagon(s) for s in [Point(ONE, ZERO), Point(*CH), Point(*CT)]]

# every eps that arises, against every eps a cut could name
seen = set()
for bits in itertools.product((0, 1), repeat=t):
    seen.add(tuple(bits[a] ^ bits[b] for a in range(t) for b in range(a + 1, t)))
allpat = set(itertools.product((0, 1), repeat=t * (t - 1) // 2))
print(f"{len(seen)} of {len(allpat)} pair-parity patterns are reachable "
      f"= 2^(t-1) = {2 ** (t - 1)}")
for e in sorted(allpat - seen):
    a, b, c = e
    assert (a + b + c) % 2 == 1, "unreachable patterns must break the triangle"
print("  every unreachable pattern breaks eps_ab + eps_bc + eps_ac = 0")

# the always-confined points, and whether their graph really is bipartite
same = {}
for a in range(t):
    for i, u in enumerate(hexes[a]):
        for j, v in enumerate(hexes[a]):
            if i < j and (i + j) % 2 == 1:
                q = Point(u.x + v.x, u.y + v.y)
                if q.norm2() not in (ONE, ZERO):
                    same[(q.x, q.y)] = (a, i, j)
pts = list(same)
adj = {i: set() for i in range(len(pts))}
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        if Point(*pts[i]).is_unit_apart(Point(*pts[j])):
            adj[i].add(j)
            adj[j].add(i)
col, ok = {}, True
for s in range(len(pts)):
    if s in col:
        continue
    col[s], stack = 0, [s]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in col:
                col[y] = 1 - col[x]
                stack.append(y)
            elif col[y] == col[x]:
                ok = False
print(f"\n{len(pts)} always-confined points at sqrt 3, "
      f"max degree {max((len(v) for v in adj.values()), default=0)}, "
      f"bipartite: {ok}")
