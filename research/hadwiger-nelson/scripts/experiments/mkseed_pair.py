"""Seed a pair search on a substrate: a = densest vertex, b = a + k e along a
unit edge at a (the one whose ray from a meets the most vertices), plus the
180-degree turn about the midpoint (swaps a and b) and both unit circles
filled along every direction the substrate uses.
usage: mkseed_pair.py <substrate.json> <out.json> <k>
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
ROOT = HN_DIR
d = json.load(open(f"{ROOT}/data/{sys.argv[1]}")); k = int(sys.argv[3])
F = Field(tuple(d["field_generators"]))
S = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(S)
key = lambda p: (round(p.fx, 9), round(p.fy, 9)); vk = {key(p): i for i, p in enumerate(S)}
ia = max(range(g.n), key=lambda v: len(g.adj[v])); a = S[ia]
def ray(u):   # how many of a + u, a + 2u, ..., a + k u are already vertices
    return sum(1 for j in range(1, k + 1) if key(a + Point(u.x * j, u.y * j)) in vk)
cands = [S[w] - a for w in g.adj[ia]]
e = max(cands, key=ray)
b = a + Point(e.x * k, e.y * k)
half = F.rational(Fr(1, 2)); m = Point((a.x + b.x) * half, (a.y + b.y) * half)
U = {}
for x, y in g.edges():
    for u in (S[y] - S[x], S[x] - S[y]): U.setdefault(key(u), u)
pts = {key(p): p for p in S}
for p in S: q = Point(m.x * 2 - p.x, m.y * 2 - p.y); pts.setdefault(key(q), q)
for u in U.values():
    pts.setdefault(key(a + u), a + u); pts.setdefault(key(b + u), b + u)
pts.setdefault(key(b), b)
P = list(pts.values())
iA = next(i for i, p in enumerate(P) if key(p) == key(a)); iB = next(i for i, p in enumerate(P) if key(p) == key(b))
G = build_graph(P)
print(f"seed {sys.argv[2]}: n={G.n} m={G.m}, pair at k={k} along a ray with {ray(e)} vertices; deg a={len(G.adj[iA])} deg b={len(G.adj[iB])}", flush=True)
json.dump({"field_generators": list(F.gens), "A": iA, "B": iB, "note": f"{sys.argv[1]}, pair a, a+{k}e, turn about the midpoint, circles filled",
           "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in P]},
          open(f"{ROOT}/data/{sys.argv[2]}", "w"))
