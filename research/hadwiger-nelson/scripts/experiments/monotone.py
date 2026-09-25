"""rho only goes DOWN when structure is added, and that reverses the strategy.

If G is an induced subgraph of W on a subset of its vertices, every proper
k-colouring of W restricts to one of G -- so a set forcing in G is forcing in
W, and

    rho(W) <= rho(G)   whenever G is contained in W.

Which means the way to a small rho is a BIGGER graph, not a smaller one. That
is the opposite of what this package assumed while the cross-pair theorem was
thought to apply to unions of copies of de Grey's G: it needs W - a - b to be
(k-1)-colourable for every cross pair, which for copies of G means G - a is
4-colourable -- exactly the vertex-criticality now retracted. The theorem is
correct; its hypothesis is not met, so unions are back on the table, and
monotonicity says they can only help.

Checked on a case small enough to compute rho exactly by brute force.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

FLD = Field((3, 11))


def moser_spindle():
    h = FLD.rational(Fraction(1, 2))
    u = Point(h, FLD.sqrt(3) * h)
    one = Point(FLD.one(), FLD.zero())
    rh = [Point(FLD.zero(), FLD.zero()), one, u, one + u]
    rot = Rotation(FLD.rational(Fraction(5, 6)),
                   FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
    return rh + [rot(p) for p in rh[1:]]


def rho_exact(g, k, cap=12):
    """Least forcing set, by brute force over subsets."""
    nv = g.n * k

    def x(v, c):
        return 1 + v * k + c

    base = [[x(v, c) for c in range(k)] for v in range(g.n)]
    for a, b in g.edges():
        for c in range(k):
            base.append([-x(a, c), -x(b, c)])

    def forcing(S):
        for c in range(k):
            with Solver(name="g4",
                        bootstrap_with=base + [[-x(v, c)] for v in S]) as s:
                if s.solve():
                    return False
        return True

    for size in range(1, min(cap, g.n) + 1):
        for S in itertools.combinations(range(g.n), size):
            if forcing(S):
                return size, S
    return None, None


sp = build_graph(moser_spindle())
r, S = rho_exact(sp, 4, cap=8)
print(f"Moser spindle, k=4: n={sp.n}, rho={r}  {S}", flush=True)

# now embed it in something larger: the spindle plus a rotated copy sharing
# the pivot, which contains the spindle as an induced subgraph
rot = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))
pts = moser_spindle()
big = list(pts)
for p in pts:
    q = rot(p)
    if q not in big:
        big.append(q)
w = build_graph(big)
print(f"  embedded in a {w.n}-vertex graph ({w.m} edges)", flush=True)
r2, S2 = rho_exact(w, 4, cap=8)
print(f"  rho(W,4) = {r2}  {S2}", flush=True)
print(f"  monotone: rho(W) <= rho(G)?  {r2} <= {r}: {r2 <= r}", flush=True)
