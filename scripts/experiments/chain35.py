"""A rhombus chain showing chi(Q(sqrt3, sqrt5)^2) >= 4.

tau = (2 + i sqrt5)/3 is a unit vector, so 4/3 = tau + conj(tau), 1/3 = tau + conj(tau) - 1 and
(1 + w)/3 = sum of the six unit vectors tau, conj(tau), -1, w tau, w conj(tau), -w (w = e^{i pi/3}).
Its squared length is |1 + w|^2/9 = 1/3. A unit rhombus of diagonal sqrt3 * u forces its tips alike
in every 3-colouring; chaining the six rhombi joins 0 to sqrt3 (1 + w)/3, a unit vector.
"""
import os, sys, json
from fractions import Fraction as Fr
HN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, HN_DIR)
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
F = Field((3, 5)); r = F.rational; s3 = F.sqrt(3); s5 = F.sqrt(5)
P = lambda x, y: Point(x, y)
tau = P(r(Fr(2, 3)), s5 * r(Fr(1, 3))); taub = P(tau.x, -tau.y)
w = P(r(Fr(1, 2)), s3 * r(Fr(1, 2)))
def mul(a, b): return P(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x)
one = P(r(1), r(0)); neg = lambda a: P(-a.x, -a.y)
us = [tau, taub, neg(one), mul(w, tau), mul(w, taub), neg(w)]
for u in us: assert u.x * u.x + u.y * u.y == F.one()
pts = [P(r(0), r(0))]
cur = pts[0]
for u in us:
    nxt = P(cur.x + s3 * u.x, cur.y + s3 * u.y)
    mid = P((cur.x + nxt.x) * r(Fr(1, 2)), (cur.y + nxt.y) * r(Fr(1, 2)))
    perp = P(-u.y * r(Fr(1, 2)), u.x * r(Fr(1, 2)))
    pts += [P(mid.x + perp.x, mid.y + perp.y), P(mid.x - perp.x, mid.y - perp.y), nxt]
    cur = nxt
end = cur
print("end point squared length:", end.x * end.x + end.y * end.y)
V = list(dict.fromkeys(pts)); g = build_graph(V); E = sorted(g.edges())
print(len(V), "vertices,", len(E), "unit-distance pairs (exact)")
def colourable(n, E, k):
    X = lambda v, c: 1 + v * k + c
    s = Solver(name="cd19")
    for v in range(n): s.add_clause([X(v, c) for c in range(k)])
    for a, b in E:
        for c in range(k): s.add_clause([-X(a, c), -X(b, c)])
    res = s.solve(); s.delete(); return res
print("3-colourable:", colourable(len(V), E, 3), " 4-colourable:", colourable(len(V), E, 4))
enc = lambda e: [[c.numerator, c.denominator] for c in e.c]
json.dump({"field_generators": [3, 5], "points": [[enc(p.x), enc(p.y)] for p in V],
           "note": "six unit rhombi of diagonal sqrt3*u for u = tau, conj(tau), -1, w tau, w conj(tau), -w; the chain joins 0 to the unit vector sqrt3(1+w)/3"},
          open(os.path.join(HN_DIR, "data", "chain35.json"), "w"))
