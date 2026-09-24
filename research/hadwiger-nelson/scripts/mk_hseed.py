"""five_rho7 u H (E-I's H turned onto the module, A at the densest vertex), with H's 446 two-edges
listed as targets.  The forced-pair search on (A, B) then only has to make those two-edges hold
'enough': H with its two-edges forces c(A) = c(B) (re-verified with three solvers)."""
import json, sys
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d7 = json.load(open(f"{ROOT}/data/five_rho7.json")); dK = json.load(open(f"{ROOT}/data/K_rot.json"))
F = Field(tuple(d7["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
P7 = [mk(xy) for xy in d7["points"]]; KP = [mk(xy) for xy in dK["points"]]
H = KP[:214]; A, B = dK["A"], dK["B"]
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
out = []; idx = {}
for p in P7 + H:
    if key(p) not in idx: idx[key(p)] = len(out); out.append(p)
Hi = [idx[key(p)] for p in H]
two = [(Hi[i], Hi[j]) for i, j in dK["two_edges"] if i < 214 and j < 214]
ser = lambda z: [[[t.numerator, t.denominator] for t in z.x.c], [[t.numerator, t.denominator] for t in z.y.c]]
json.dump({"field_generators": list(F.gens), "A": Hi[A], "B": Hi[B], "two_edges": two,
           "note": "five_rho7 u H (E-I, turned 150 degrees, A at the densest vertex); two_edges = H's 446 distance-2 pairs",
           "points": [ser(z) for z in out]}, open(f"{ROOT}/data/rho7_H.json", "w"))
print(len(out), "points;", len(two), "two-edges; H points already in five_rho7:", sum(1 for p in H if idx[key(p)] < len(P7)))
