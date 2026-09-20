"""What is special about the five points that meet every colour class?

Sa admits at least forty thousand distinct realisable colour classes, each
about a quarter of the graph.  A random five points would miss thousands of
them.  The witness [49, 91, 210, 215, 330] misses none.

So the five are not scattered -- they are placed.  This reports what they look
like: their degrees, their mutual distances, what they induce, and how they
sit relative to the pivot structure Sa is built around.
"""
import sys, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.graph import build_graph
from pysat.solvers import Solver

g = build_graph(degrey.build_Sa())
S = [49, 91, 210, 215, 330]
deg = g.degrees()
V = g.vertices
print(f"Sa: {g.n} vertices; witness {S}", flush=True)
print(f"  degrees {[deg[v] for v in S]} against max {max(deg)}", flush=True)
sub = g.induced(sorted(S))
print(f"  induced: {sub.n} vertices, {sub.m} edges", flush=True)
print("  pairwise squared distances:")
for a, b in itertools.combinations(S, 2):
    d = V[a].dist2(V[b])
    print(f"    {a:4} {b:4}  {float(d):8.4f}"
          + ("   unit" if float(d) == 1.0 else ""), flush=True)
print(f"  distances from the origin: "
      f"{[round(float(V[v].norm2()) ** 0.5, 3) for v in S]}", flush=True)
hub = max(range(g.n), key=lambda v: deg[v])
print(f"  distances from the hub (vertex {hub}, degree {deg[hub]}): "
      f"{[round(float(V[v].dist2(V[hub])) ** 0.5, 3) for v in S]}", flush=True)
print(f"  how many share a neighbour with the hub: "
      f"{sum(1 for v in S if g.adj[v] & g.adj[hub])}", flush=True)

# how much does each one matter?  drop it and see how many classes escape
k = 4


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
print("\n  dropping each point, is the rest still forcing?")
for drop in S:
    T = [v for v in S if v != drop]
    esc = None
    for c in range(k):
        with Solver(name="cd19",
                    bootstrap_with=cls + [[-x(v, c)] for v in T]) as s:
            if s.solve():
                esc = c
                break
    verdict = "forcing still" if esc is None else f"colour {esc} escapes"
    print(f"    without {drop}: {verdict}", flush=True)
