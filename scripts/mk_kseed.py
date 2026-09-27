"""Exoo-Ismailescu's K inside our modules.

E-I's 30 unit vectors, turned by 150 degrees, are unit vectors of five_247_c (hence of five_rho7
and of the blocked module).  So their graph H, turned by 150 degrees and put at a vertex O of
five_rho7, lies in five_rho7's module, and K = H u lambda_A(H) lies in five_rho7 u lambda_O(five_rho7).
Writes: data/K_rot.json (K alone, with its 892 two-edges listed) and
        data/rho7_lambda_K.json (five_rho7 u lambda_O(five_rho7) u K, A = O)."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import json, sys
from fractions import Fraction as Fr
sys.path.insert(0, HN_DIR)
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
ROOT = HN_DIR
d7 = json.load(open(f"{ROOT}/data/five_rho7.json")); F = Field(tuple(d7["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
P7 = [mk(xy) for xy in d7["points"]]
dH = json.load(open(f"{ROOT}/data/ei_H214.json")); FH = Field(tuple(dH["field_generators"]))
# lift H into F (same basis names: F contains sqrt3, sqrt11)
basisF = [F.element([Fr(1) if j == i else Fr(0) for j in range(len(F.one().c))]) for i in range(len(F.one().c))]
basisH = [FH.element([Fr(1) if j == i else Fr(0) for j in range(len(FH.one().c))]) for i in range(len(FH.one().c))]
# identify FH basis elements with F elements by their squares (1, sqrt3, sqrt11, sqrt33)
def lift_basis(eH):
    sq = (eH * eH).c[0]
    for eF in basisF:
        if eF * eF == F.rational(sq): return eF
    raise ValueError
imgs = [lift_basis(e) for e in basisH]
lift = lambda z: sum((imgs[i] * F.rational(c) for i, c in enumerate(z.c)), F.zero())
H = [Point(lift(FH.element([Fr(a, b) for a, b in xy[0]])), lift(FH.element([Fr(a, b) for a, b in xy[1]]))) for xy in dH["points"]]
A, B = dH["A"], dH["B"]
r3 = next(e for e in basisF if e * e == F.rational(3)); r11 = next(e for e in basisF if e * e == F.rational(11))
q = lambda a, b: F.rational(Fr(a, b))
R150 = Rotation(-r3 * q(1, 2), q(1, 2))                         # cos 150 = -sqrt3/2, sin 150 = 1/2
lam = Rotation(q(49, 50), r11 * q(3, 50))                       # E-I's lambda: cos 49/50, sin 3 sqrt11/50
g7 = build_graph(P7)
deg = [0] * len(P7)
for a, b in g7.edges(): deg[a] += 1; deg[b] += 1
O = max(range(len(P7)), key=lambda i: deg[i])
Hr = [P7[O] + R150(p - H[A]) for p in H]                       # A goes to O
lamO = lam.about(P7[O])
K = Hr + [lamO(p) for p in Hr]
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
Ku = []; seen = {}
for p in K:
    if key(p) not in seen: seen[key(p)] = len(Ku); Ku.append(p)
iA = seen[key(Hr[A])]; iB = seen[key(Hr[B])]; iB2 = seen[key(lamO(Hr[B]))]
two = [(i, j) for i in range(len(Ku)) for j in range(i + 1, len(Ku)) if Ku[i].dist2(Ku[j]) == 4]
print(f"K: {len(Ku)} points; A={iA}, B={iB}, B'={iB2}, |BB'|^2 = {Ku[iB].dist2(Ku[iB2])}; two-edges {len(two)}")
ser = lambda z: [[[t.numerator, t.denominator] for t in z.x.c], [[t.numerator, t.denominator] for t in z.y.c]]
json.dump({"field_generators": list(F.gens), "A": iA, "B": iB, "Bp": iB2, "two_edges": two,
           "note": "E-I K = H u lambda_A(H), H turned 150 degrees onto five_rho7's units, A at five_rho7's densest vertex",
           "points": [ser(z) for z in Ku]}, open(f"{ROOT}/data/K_rot.json", "w"))
# check H's unit edges are unit edges of five_rho7's unit set
U7 = {key(P7[b] - P7[a]) for a, b in g7.edges()} | {key(P7[a] - P7[b]) for a, b in g7.edges()}
gH = build_graph(Hr)
print("H turned: its unit edges along five_rho7 units:", all(key(Hr[b] - Hr[a]) in U7 for a, b in gH.edges()), f"({gH.m} edges)")
full = list(P7) + [lamO(p) for p in P7] + Ku
out = []; s2 = {}
for p in full:
    if key(p) not in s2: s2[key(p)] = len(out); out.append(p)
json.dump({"field_generators": list(F.gens), "A": s2[key(Hr[A])], "B": s2[key(Hr[B])], "Bp": s2[key(lamO(Hr[B]))],
           "two_edges": [(s2[key(Ku[i])], s2[key(Ku[j])]) for i, j in two],
           "note": "five_rho7 u lambda_O(five_rho7) u K (E-I's graph, turned onto the module), O = A",
           "points": [ser(z) for z in out]}, open(f"{ROOT}/data/rho7_lambda_K.json", "w"))
print(f"rho7_lambda_K: {len(out)} points")
