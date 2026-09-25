"""Kempe component size as a rigidity measure -- cheap, and never used here.

A Kempe swap exchanges two colours on one connected component of the subgraph
they induce, and the components' sizes say how FREE a colouring is: in a loose
graph the two-colour subgraphs shatter into small pieces and a colouring can be
edited locally; in a rigid one they percolate, and changing one vertex's colour
drags a quarter of the graph with it.

On the tight hexagon graph the random walk averaged about a thousand vertices
per swapped component, out of 4159 -- a quarter of the graph.  That number means
nothing without a baseline, so compute the same statistic, from one colouring
and every (vertex, colour) pair, on the graphs it was built from.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import deque
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
for NAME in sys.argv[1:]:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    G = build_graph(P); n = G.n; E = list(G.edges())
    adj = [[] for _ in range(n)]
    for a, b in E: adj[a].append(b); adj[b].append(a)
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K): cnf.append([-X(v, a), -X(v, b)])
    for a, b in E:
        for c in range(K): cnf.append([-X(a, c), -X(b, c)])
    s = Solver(name="cd19", bootstrap_with=cnf); ok = s.solve()
    m = s.get_model(); s.delete()
    col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
    # every Kempe component, for every pair of colours
    sizes = []
    for a in range(K):
        for b in range(a + 1, K):
            seen = set()
            for v in range(n):
                if col[v] not in (a, b) or v in seen: continue
                comp = 0; q = deque([v]); seen.add(v)
                while q:
                    u = q.popleft(); comp += 1
                    for w in adj[u]:
                        if w not in seen and col[w] in (a, b):
                            seen.add(w); q.append(w)
                sizes.append(comp)
    # the size a uniformly random vertex's component has (size-biased mean)
    tot = sum(sizes); sb = sum(x * x for x in sizes) / tot
    big = max(sizes)
    print(f"  {NAME:<28s} n={n:<5d} components {len(sizes):<5d} "
          f"largest {big:<5d} ({100*big/n:.0f}% of n)  "
          f"size-biased mean {sb:7.1f} ({100*sb/n:.1f}% of n)   "
          f"[{time.time()-t0:.0f}s]", flush=True)
