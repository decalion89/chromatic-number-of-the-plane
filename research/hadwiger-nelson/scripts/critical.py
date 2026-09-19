"""Is de Grey's G really 5-vertex-critical?  The claim was asserted, not tested.

A k-vertex-critical graph has chi(G - v) < k for EVERY v.  This package leans
on G being 5-vertex-critical to explain why nothing is ever forced at five
colours, and the evidence offered was separability -- which is a CONSEQUENCE
of criticality, not a proof of it.

There is a strong reason to doubt it: Heule reduced de Grey's graph to 874 and
then to 553 vertices, so 5-chromatic proper subgraphs exist.  If one of them
sits inside G, then G minus a vertex outside it is still 5-chromatic and G is
not vertex-critical at all.

The cheap decisive test is a proper SUBGRAPH: if any k-core of G fails to be
4-colourable, G has a 5-chromatic proper subgraph and criticality is dead.
Cores come free, and the SAT call is the same shape as the one that certified
the graph in the first place.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
deg = g.degrees()
print(f"G: {g.n} vertices, {g.m} edges, degrees "
      f"{min(deg)}..{max(deg)}", flush=True)

for kc in (4, 5, 6, 7, 8, 9, 10):
    core = g.k_core(kc)
    print(f"  {kc}-core: {core.n} vertices, {core.m} edges", flush=True)
    if core.n == 0:
        break
    if core.n == g.n:
        continue
    # a proper subgraph: is it still non-4-colourable?
    k = 4
    cls = [[1 + v * k + c for c in range(k)] for v in range(core.n)]
    for a, b in core.edges():
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    t = time.time()
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"    4-colourable: {r}  [{time.time()-t:.0f}s]", flush=True)
    if not r:
        print(f"    *** a PROPER subgraph of {core.n} vertices is already "
              f"5-chromatic: G is NOT 5-vertex-critical ***", flush=True)
        break
