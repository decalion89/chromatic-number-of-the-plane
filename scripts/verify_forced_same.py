"""Verify a claimed forced pair c(A) = c(B) at |AB| = 5 and its lambda-closure, from scratch.

1. rebuild the graph from exact coordinates with build_graph (ALL unit distances, exact);
2. re-check every edge as an exact identity |p - q|^2 = 1, and |AB|^2 = 25;
3. 'c(A) != c(B)' must be UNSAT for cadical AND glucose (no triangle pinned);
4. close: lambda = (49 + 3 sqrt-11)/50 about A sends B to B' with |BB'|^2 = 1 (exact);
   the union G u lambda_A(G), rebuilt exactly, must be UNSAT for plain 5-colouring
   with two solvers -- that union is the 6-chromatic unit-distance graph.
Nothing is claimed unless every step agrees."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, HN_DIR)
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver
t0 = time.time(); K = 5
d = json.load(open(sys.argv[1])); F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; A, B = d["A"], d["B"]
g = build_graph(V); E = list(g.edges()); one = F.rational(1)
assert all((V[a] - V[b]).norm2() == one for a, b in E)
print(f"G: n={g.n} m={len(E)} (exact); |AB|^2 = {V[A].dist2(V[B])}   [{time.time()-t0:.0f}s]", flush=True)
assert V[A].dist2(V[B]) == F.rational(25)
X = lambda v, c: 1 + v * K + c
def cnf(n, E):
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K): cl.append([-X(a, c), -X(b, c)])
    return cl
base = cnf(g.n, E) + [[-X(A, c), -X(B, c)] for c in range(K)]
for name in ("cadical195", "glucose4"):
    s = Solver(name=name, bootstrap_with=base); r = s.solve(); s.delete()
    print(f"  {name}: c(A) != c(B) satisfiable? {r}   [{time.time()-t0:.0f}s]", flush=True)
    if r: print("NOT FORCED -- stop"); sys.exit(1)
basis = [F.element([Fr(1) if j == i else Fr(0) for j in range(len(F.one().c))]) for i in range(len(F.one().c))]
r11 = next(e for e in basis if e * e == F.rational(11))
lam = Rotation(F.rational(Fr(49, 50)), r11 * F.rational(Fr(3, 50))).about(V[A])
W = [lam(p) for p in V]
print(f"  |B lambda(B)|^2 = {V[B].dist2(W[B])}", flush=True)
assert V[B].dist2(W[B]) == one
U = list({p: p for p in V + W}.values())
gu = build_graph(U); Eu = list(gu.edges())
assert all((U[a] - U[b]).norm2() == one for a, b in Eu)
print(f"union: n={gu.n} m={len(Eu)} (exact)   [{time.time()-t0:.0f}s]", flush=True)
for name in ("cadical195", "glucose4"):
    s = Solver(name=name, bootstrap_with=cnf(gu.n, Eu)); r = s.solve(); s.delete()
    print(f"  {name}: union 5-colourable? {r}   [{time.time()-t0:.0f}s]", flush=True)
