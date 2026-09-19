"""What the confined set looks like, since its degeneracy will not move.

Maximum degree 4 and degeneracy 2 is a specific shape, not a generic one: it
says every subgraph has a vertex of degree at most two, so the graph peels
down to nothing two edges at a time.  Unions of paths and cycles do that, and
so do cacti, and none of them can be 4-chromatic.  This reports the degree
distribution, the components and their cycle structure on de Grey's own
configuration, to say which it is.
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
CH = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
              FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))

T = 3
seeds, q = [], Point(ONE, ZERO)
for _ in range(T):
    seeds.append(q)
    q = CH(q)
hexes = []
for s in seeds:
    h, p = [], s
    for _ in range(6):
        h.append(p)
        p = ROT60(p)
    hexes.append(h)

aux = {}
for a in range(T):
    for b in range(a, T):
        for i, u in enumerate(hexes[a]):
            for j, v in enumerate(hexes[b]):
                if (a, i) >= (b, j):
                    continue
                s = Point(u.x + v.x, u.y + v.y)
                if s.norm2() in (ONE, ZERO):
                    continue
                aux.setdefault((s.x, s.y), set()).add((a, i))
                aux[(s.x, s.y)].add((b, j))
keys = list(aux)
adj = {i: set() for i in range(len(keys))}
for i in range(len(keys)):
    pi = Point(*keys[i])
    for j in range(i + 1, len(keys)):
        if pi.is_unit_apart(Point(*keys[j])):
            adj[i].add(j)
            adj[j].add(i)

for bits in itertools.product((0, 1), repeat=T):
    col = {(a, i): (i + bits[a]) % 2 for a in range(T) for i in range(6)}
    conf = {i for i in range(len(keys)) if {col[n] for n in aux[keys[i]]} == {0, 1}}
    sub = {v: adj[v] & conf for v in conf}
    m = sum(len(s) for s in sub.values()) // 2
    dist = {}
    for v in conf:
        dist[len(sub[v])] = dist.get(len(sub[v]), 0) + 1
    # components, and the cycle rank m - n + c of each
    seen, comps = set(), []
    for v in conf:
        if v in seen:
            continue
        stack, part = [v], []
        seen.add(v)
        while stack:
            x = stack.pop()
            part.append(x)
            for y in sub[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        e = sum(len(sub[x]) for x in part) // 2
        comps.append((len(part), e, e - len(part) + 1))
    ranks = sorted({r for _, _, r in comps})
    print(f"bits {bits}: {len(conf)} confined, {m} edges, "
          f"degrees {dict(sorted(dist.items()))}, {len(comps)} components, "
          f"cycle ranks {ranks}", flush=True)
