"""The 803-graph with rho7-rotated copies about its densest vertex.

rho7 = (1 + 4 sqrt-3)/7 is the turn between two Eisenstein vectors of norm 7:
cos 1/7, sin 4 sqrt3/7, rational in Q(sqrt3), integral at 5.  It keeps the
field and the 803-module's divisibility by 2, 3, 4 and 8 (the blocking a
fifth colour needs), and -- measured on small modules -- is exactly what kills
the twisted coset colourings.
"""
import sys, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
c = P[max(range(g.n), key=lambda v: len(g.adj[v]))]
r3 = F.sqrt(3)
rots = {"rho7": (F.rational(Fr(1, 7)), r3 * F.rational(Fr(4, 7))), "rho7inv": (F.rational(Fr(1, 7)), -r3 * F.rational(Fr(4, 7)))}
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
pts = {key(p): p for p in P}
for name in sys.argv[2:]:
    ca, sa = rots[name]
    for p in P:
        dx, dy = p.x - c.x, p.y - c.y
        q = Point(c.x + dx * ca - dy * sa, c.y + dx * sa + dy * ca); pts.setdefault(key(q), q)
V = list(pts.values()); G = build_graph(V)
print(f"803 + {sys.argv[2:]}: n={G.n} m={G.m}", flush=True)
json.dump({"field_generators": list(F.gens), "note": f"five_247_c u {sys.argv[2:]} about its densest vertex",
           "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in V]},
          open(f"{ROOT}/data/{sys.argv[1]}", "w"))
