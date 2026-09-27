"""Write the blocked substrate: the 803-graph and its lambda-rotated copy.

lambda = (49 + 3 sqrt-11)/50 is Exoo-Ismailescu's rotation (the 5,5,1 triangle);
its denominator carries 5, which makes the edge module 5-divisible at some edge
vectors and kills every coset colouring mod 5 -- measured: blocks at 2, 3, 4, 5.
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
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
v0 = max(range(g.n), key=lambda v: len(g.adj[v])); c = P[v0]
ca = F.rational(Fr(49, 50)); sa = F.rational(Fr(3, 50)) * F.sqrt(11)
def rot(p):
    dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + dx * ca - dy * sa, c.y + dx * sa + dy * ca)
U = list(dict.fromkeys(P + [rot(p) for p in P]))
G = build_graph(U)
json.dump({"field_generators": list(F.gens), "note": "five_247_c u lambda_c(five_247_c), lambda=(49+3sqrt-11)/50 about the densest vertex",
           "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in U]},
          open(f"{ROOT}/data/five_247_c_lambda.json", "w"))
print("written", G.n, G.m)
