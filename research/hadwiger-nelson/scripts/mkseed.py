"""Seed a blocked gadget search: substrate u rotated copies about its densest
vertex (rotations with 5 in the denominator), pair a = that vertex and
b = a + 2e along a rotated edge, plus the 180-degree turn about the midpoint and
both unit circles filled along every direction.

usage: mkseed.py <substrate.json> <out.json> <rotation names...>   (lam, nu, nubar)
"""
import sys, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/{sys.argv[1]}"))
F = Field(tuple(d["field_generators"]))
G0 = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(G0)
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
c0 = max(range(g.n), key=lambda v: len(g.adj[v])); c = G0[c0]
r11 = F.sqrt(11)
ROT = {"lam": (F.rational(Fr(49, 50)), F.rational(Fr(3, 50)) * r11),
       "nu": (F.rational(Fr(-1, 10)), F.rational(Fr(3, 10)) * r11),
       "nubar": (F.rational(Fr(-1, 10)), -F.rational(Fr(3, 10)) * r11)}
def rot(p, cs):
    ca, sa = cs; dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + dx * ca - dy * sa, c.y + dx * sa + dy * ca)
pts = {key(p): p for p in G0}
first = None
for name in sys.argv[3:]:
    cs = ROT[name]
    for p in G0: pts.setdefault(key(rot(p, cs)), rot(p, cs))
    if first is None: first = cs
S = list(pts.values())
gs = build_graph(S)
ia = next(i for i, p in enumerate(S) if key(p) == key(c))
# pair direction: an edge at a that came from the first rotated copy (5 in its denominator)
cands = [S[w] - c for w in gs.adj[ia]]
e = max(cands, key=lambda u: max(t.denominator for t in list(u.x.c) + list(u.y.c)))
a, b = c, c + e + e; m = c + e
U = {}
for x, y in gs.edges():
    for u in (S[y] - S[x], S[x] - S[y]): U.setdefault(key(u), u)
turn = [Point(m.x * 2 - p.x, m.y * 2 - p.y) for p in S]
circ = [a + u for u in U.values()] + [b + u for u in U.values()]
for p in turn + circ + [b]: pts.setdefault(key(p), p)
P = list(pts.values())
ia = next(i for i, p in enumerate(P) if key(p) == key(a)); ib = next(i for i, p in enumerate(P) if key(p) == key(b))
G = build_graph(P)
print(f"seed: n={G.n} m={G.m}, {len(U)} directions, deg a={len(G.adj[ia])} deg b={len(G.adj[ib])}, e={e}", flush=True)
json.dump({"field_generators": list(F.gens), "A": ia, "B": ib, "note": f"{sys.argv[1]} + {sys.argv[3:]} about its densest vertex",
           "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in P]},
          open(f"{ROOT}/data/{sys.argv[2]}", "w"))
