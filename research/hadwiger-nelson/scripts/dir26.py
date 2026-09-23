import sys, json, math
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
E = edge_vectors(g); B = echelon(E); C = [coords(B, v) for v in E]
raw = {}
for a, b in g.edges():
    dx = g.vertices[b].x - g.vertices[a].x; dy = g.vertices[b].y - g.vertices[a].y
    raw[tuple(dx.c) + tuple(dy.c)] = (dx, dy)
den = 1
for v in raw:
    for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
for k in (26,):
    v = E[k]
    for key, (dx, dy) in raw.items():
        iv = tuple(int(Fr(q) * den) for q in key)
        if iv == v or tuple(-x for x in iv) == v:
            print(f"direction {k}: e = ({dx}, {dy})  ~ ({float(dx):.6f}, {float(dy):.6f}), angle {math.degrees(math.atan2(float(dy), float(dx))):.4f} deg")
            break
    print("  divisibility in M:", reduce(gcd, [abs(t) for t in C[k]]))
    cnt = sum(1 for a, b in g.edges() if tuple(int(Fr(q) * den) for q in tuple((g.vertices[b].x - g.vertices[a].x).c) + tuple((g.vertices[b].y - g.vertices[a].y).c)) in (v, tuple(-x for x in v)))
    print("  edges of the 803-graph along it:", cnt)
