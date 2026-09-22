"""The ladder: cap -> forced pair -> contradiction.

G is Ya u Yb.  The measurement says those two 791-point halves share exactly
one vertex and are joined by exactly one edge.  A single edge cannot raise the
chromatic number of a disjoint union by itself, and neither can a shared
vertex.  So G being 5-chromatic while Y is 4-chromatic forces a very specific
statement about Y:

    in EVERY 4-colouring of Y, the colour of the cross-edge endpoint p is
    determined by the colour of the shared vertex s.

Then in G the two halves agree at s, hence agree at p_a and p_b -- and p_a p_b
is an edge.  That is the whole contradiction.  It is the Moser spindle with Y
in the role of the rhombus: the rhombus forces its two tips equal at three
colours, Y forces (s, p) equal at four.

This script finds s and p and tests that claim outright.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Y, _rot_half_pi_pm
from hn.geometry import Point
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 5, 7, 11))
PIV = Point(F.rational(-2), F.zero())
Y = build_Y(F)
Ya = [_rot_half_pi_pm(F, +1).about(PIV)(p) for p in Y]
Yb = [_rot_half_pi_pm(F, -1).about(PIV)(p) for p in Y]
sa, sb = set(Ya), set(Yb)
shared = sa & sb
print(f"Y: {len(Y)} points;  Ya n Yb = {len(shared)}", flush=True)
s = next(iter(shared))
print(f"  shared vertex s = {s.approx()}", flush=True)
cross = [(i, j) for i, p in enumerate(Ya) if p not in shared
         for j, q in enumerate(Yb) if q not in shared and p.is_unit_apart(q)]
print(f"  cross edges: {len(cross)}", flush=True)
for i, j in cross:
    print(f"    Ya[{i}] {Ya[i].approx()}  --  Yb[{j}] {Yb[j].approx()}", flush=True)

# pull s and the cross endpoint back into Y's own coordinates
back_a = _rot_half_pi_pm(F, +1).about(PIV).__self__ if False else None
rota, rotb = _rot_half_pi_pm(F, +1).about(PIV), _rot_half_pi_pm(F, -1).about(PIV)
idx = {p: i for i, p in enumerate(Y)}
i0, j0 = cross[0]
# Ya[i] = rota(Y[i]) because the list order is preserved
s_in_Y_a = Y[Ya.index(s)]
s_in_Y_b = Y[Yb.index(s)]
p_a, p_b = Y[i0], Y[j0]
print(f"\n  in Y's own coordinates:")
print(f"    s  (via Ya) = {s_in_Y_a.approx()}      s (via Yb) = {s_in_Y_b.approx()}")
print(f"    p  (via Ya) = {p_a.approx()}           p (via Yb) = {p_b.approx()}")
print(f"    |s - p|^2 (a-side) = {(s_in_Y_a - p_a).norm2()}")
print(f"    |s - p|^2 (b-side) = {(s_in_Y_b - p_b).norm2()}", flush=True)

g = build_graph(Y)
n = g.n
pos = {p: i for i, p in enumerate(g.vertices)}
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])

def forced_equal(u, v):
    """Can u and v differ?  By colour permutation it is enough to test
    c(u)=0, c(v)=1."""
    s2 = Solver(name="cd15", bootstrap_with=cnf)
    r = s2.solve(assumptions=[X(u, 0), X(v, 1)])
    s2.delete(); return not r

def forced_diff(u, v):
    s2 = Solver(name="cd15", bootstrap_with=cnf)
    r = s2.solve(assumptions=[X(u, 0), X(v, 0)])
    s2.delete(); return not r

for tag, su, pu in (("a-side", s_in_Y_a, p_a), ("b-side", s_in_Y_b, p_b)):
    u, v = pos[su], pos[pu]
    t0 = time.time()
    fe, fd = forced_equal(u, v), forced_diff(u, v)
    print(f"\n  {tag}: s=v{u} p=v{v}   forced EQUAL: {fe}   forced DIFFERENT: {fd}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
