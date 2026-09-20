"""Settle one vertex: is G - 1420 still 5-chromatic?

v = 1420 has degree 4 and its removal made the 4-colourability question hard
enough to exhaust a three-million-conflict budget, while every essential
vertex answered quickly with a colouring.  That asymmetry is the signature of
UNSAT, so this runs the single instance to completion with no budget.

UNSAT means G - v is still 5-chromatic, which means G is NOT 5-vertex-critical
-- and the corollary this package leans on does not apply to G.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
v = 1420
print(f"G: {g.n} vertices; removing v={v} (degree {len(g.adj[v])})", flush=True)
sub = g.induced([u for u in range(g.n) if u != v])
k = 4
cls = [[1 + i * k + c for c in range(k)] for i in range(sub.n)]
for a, b in sub.edges():
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
# the same pinning the certificate uses: a triangle takes 0, 1, 2
tri = sub.find_clique(3)
if tri:
    for i, u in enumerate(tri):
        cls.append([1 + u * k + i])
    print(f"  pinned triangle {tri} to colours 0,1,2", flush=True)
print(f"  {sub.n} vertices, {sub.m} edges, {len(cls)} clauses", flush=True)
t = time.time()
with Solver(name="cd19", bootstrap_with=cls) as s:
    r = s.solve()
print(f"  4-colourable: {r}  [{time.time()-t:.0f}s]", flush=True)
if r is False:
    print("  *** G - 1420 is STILL 5-chromatic: G is NOT 5-vertex-critical ***")
else:
    print("  G - 1420 became 4-colourable: v is essential after all, and the "
          "budget exhaustion was only slowness")
