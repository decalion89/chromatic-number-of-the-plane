"""Independent verification of a claimed distance-2 gadget, before a word is said.

A gadget is a unit-distance graph H with two vertices a, b at distance exactly 2
that no proper 5-colouring of H colours alike.  With Exoo-Ismailescu's theorem
(every 5-colouring of the plane has a monochromatic pair at distance 1 or 2)
that is a proof of chi(R^2) >= 6, so the claim is re-derived by a route that
shares nothing with the search that produced it:

  1. rebuild H from the saved exact coordinates and re-verify every edge as an
     exact field identity |p - q|^2 == 1;
  2. check |a - b|^2 == 4 exactly;
  3. MERGE a and b into one vertex (a graph operation, not an assumption) and
     ask three unrelated solvers whether the merged graph is 5-colourable, with
     no colour pinned anywhere;
  4. report any disagreement as a disagreement, never as a result.
"""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time(); K = 5
d = json.load(open(sys.argv[1]))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
a, b = d["A"], d["B"]
g = build_graph(P); n = g.n; E = list(g.edges())
one = F.rational(1)
bad = [(x, y) for x, y in E if (g.vertices[x] - g.vertices[y]).norm2() != one]
print(f"n={n} edges={len(E)}; edges failing the exact unit test: {len(bad)}; "
      f"|a-b|^2 = {g.vertices[a].dist2(g.vertices[b])}", flush=True)
assert not bad and g.vertices[a].dist2(g.vertices[b]) == F.rational(4)
# merge b into a
ren = lambda v: a if v == b else (v if v < b else v - 1)
ME = {(min(ren(x), ren(y)), max(ren(x), ren(y))) for x, y in E}
assert all(x != y for x, y in ME), "a and b adjacent?"
m = n - 1
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(m)]
for x, y in ME:
    for c in range(K): cnf.append([-X(x, c), -X(y, c)])
for name in ("cd19", "g4", "m22"):
    s = Solver(name=name, bootstrap_with=cnf)
    r = s.solve(); s.delete()
    print(f"  {name}: merged graph 5-colourable = {r}   [{time.time()-t0:.0f}s]", flush=True)
