"""Why rho(Sa, 4) is 7, and what that says about the distance left.

A k-chromatic subgraph uses all k colours in every k-colouring of the whole
graph, so it is rainbow-forcing outright:

    THEOREM.  rho(G, k) <= |H| for every k-chromatic subgraph H of G.

de Grey's Sa is 4-chromatic and contains Moser spindles, which are 4-critical
on seven vertices -- so rho(Sa, 4) <= 7, and the measurement said exactly 7.
This checks that the seven the construction found really are one.

The consequence is the honest measure of what is left.  A core of three at a
pivot of degree 60 needs rho <= 63, and the way to a small rho is a small
k-chromatic subgraph.  At four colours that is the Moser spindle, seven
vertices, and everything works.  At five the smallest published
5-chromatic unit-distance graph has around five hundred vertices.  Seven
against sixty-three is comfortable; five hundred against sixty-three is not.
"""
import sys
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.graph import build_graph
from hn.forced import ColourRelations, forcing_set, shrink_forcing_set
from pysat.solvers import Solver

g = build_graph(degrey.build_Sa())
rel = ColourRelations(g, 4)
S, closed = forcing_set(rel, limit=500)
S = shrink_forcing_set(rel, S) if closed else S
print(f"rho(Sa, 4) = {len(S)}, closed={closed}", flush=True)

sub = g.induced(sorted(S))
deg = sub.degrees()
print(f"  induced: {sub.n} vertices, {sub.m} edges, degrees {sorted(deg)}")


def chi(gr, hi=5):
    for k in range(1, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(gr.n)]
        for a, b in gr.edges():
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return None


c = chi(sub)
print(f"  chromatic number of the induced subgraph: {c}")
crit = all(chi(sub.induced([u for u in range(sub.n) if u != v])) < c
           for v in range(sub.n))
print(f"  vertex-critical: {crit}")
if sub.n == 7 and sub.m == 11 and c == 4 and crit:
    print("  => it is a Moser spindle: 7 vertices, 11 edges, 4-critical")
    print("     so rho(Sa,4) <= 7 was forced by the subgraph, and is exactly 7")
print("\nthe distance left, stated through the same theorem:")
print("  a core of three at a degree-60 pivot needs rho <= 63")
print("  at k=4 the smallest 4-chromatic unit-distance graph has 7 vertices")
print("  at k=5 the smallest published one has around 500")
