"""Does a rotation that is NOT integral at 5 escape the coset colourings?

RESIDUE_DEGREE_THREE assumed every edge vector integral at 5.  In the Moser
field 5 splits in Q(sqrt-11) into two primes swapped by complex conjugation, so
a unit vector can carry 5 in its denominator: lambda = (49 + 3 sqrt-11)/50 has
valuations +2 and -2 there -- and it is exactly Exoo-Ismailescu's rotation.
Screen: E-I's H (integral) and K = H u lambda_A(H), unit edges only; the
803-graph, and the 803-graph with a lambda-rotated copy.
"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.homcol import screen
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/ei_rebuild.py").read().split("E1, E2 = graph(V)")[0])
extra = [(-2,0,0,-6),(8,0,0,4),(-4,-6,-6,-4),(-4,6,6,-4),(-3,-3,-3,-5),(-4,0,-12,4),
         (-4,0,12,4),(7,-3,3,3),(7,3,-3,3)]
A = P(*extra[0])
VH = list(dict.fromkeys(V + [P(*e) for e in extra]))
ca, sa = F.rational(Fr(49, 50)), F.rational(Fr(3, 50)) * r11
def rot_about(c, p):
    dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + dx * ca - dy * sa, c.y + dx * sa + dy * ca)
VK = list(dict.fromkeys(VH + [rot_about(A, p) for p in VH]))
for name, pts in (("E-I H, unit edges", VH), ("E-I K = H u lambda(H), unit edges", VK)):
    g = build_graph(pts); r = screen(g, 5)
    print(f"  {name}: n={g.n} m={g.m}, {r['edge_vectors']} edge vectors dim {r['dimension']}: {r['verdict']}", flush=True)
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F2 = Field(tuple(d["field_generators"]))
P2 = [Point(F2.element([Fr(a, b) for a, b in x]), F2.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P2)
v0 = max(range(g.n), key=lambda v: len(g.adj[v])); c = P2[v0]
ca = F2.rational(Fr(49, 50)); sa = F2.rational(Fr(3, 50)) * F2.sqrt(11)
g2 = build_graph(list(dict.fromkeys(P2 + [rot_about(c, p) for p in P2])))
r = screen(g2, 5)
print(f"  803 u lambda-rotated copy: n={g2.n} m={g2.m}, {r['edge_vectors']} edge vectors dim {r['dimension']}: {r['verdict']}", flush=True)
