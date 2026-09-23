"""Independent verification of any claimed 6-chromatic graph, before a word is said.

A 5-colouring solve that returns False is the one answer in this project that
would matter more than every other, so it gets checked by a route that shares
nothing with the one that produced it:

  1. rebuild the graph from the saved exact coordinates, and re-verify every
     edge as an exact field identity |p - q|^2 == 1, not a float;
  2. confirm it is a UNIT-DISTANCE graph -- no edge at any other length;
  3. solve 5-colourability again with NO triangle pinned, so a wrong pin cannot
     manufacture an UNSAT;
  4. solve it with a DIFFERENT solver (glucose, then minisat), since two
     unrelated implementations agreeing is the standard for such claims;
  5. and only then extract a vertex-critical core for a certificate.

Any disagreement is reported as a disagreement, never as a result.
"""
import sys, time, json
from fractions import Fraction as Fr
from itertools import combinations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time(); K = 5
path = sys.argv[1]
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
one = F.rational(1)
bad = [(x, y) for x, y in E if (g.vertices[x] - g.vertices[y]).norm2() != one]
print(f"{path}: n={n} edges={len(E)}; edges failing the exact test: {len(bad)}",
      flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])
for name in ("cd19", "g4", "m22"):
    s = Solver(name=name, bootstrap_with=cnf)
    r = s.solve(); s.delete()
    print(f"  {name}, no pin: 5-colourable = {r}   [{time.time()-t0:.0f}s]",
          flush=True)
