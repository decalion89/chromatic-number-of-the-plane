"""Could a core of three force STRUCTURALLY?  Only by breaking a record.

A core of three at p needs N(p) u T forcing, and there are two ways to force:
contain a 5-chromatic subgraph, or be forced by the ambient graph.  Take the
first.

    If N(p) u T contains a 5-chromatic subgraph, then N(p) u T is itself a
    5-chromatic unit-distance graph on deg(p) + 3 vertices -- 63 at de Grey's
    best pivot, against a published record near five hundred.

So the structural route to a core of three is not merely hard: it would smash
the smallest-5-chromatic-graph record by a factor of eight as a side effect.
It is worth checking how close the geometry even gets.

The circle is bipartite -- its own unit-distance graph is paths and 6-cycles --
and each target is one away from at most TWO circle points, because two unit
circles meet twice.  So colour the circle with 2 and each target sees at most
two colours; three targets forming a triangle would need three more, giving 5.
That is the only shape that could work, and this looks for it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.solvers import Solver
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph


def chi(g, verts, hi=5):
    vs = sorted(verts)
    pos = {v: i for i, v in enumerate(vs)}
    es = [(pos[a], pos[b]) for a in vs for b in g.adj[a] if b in pos and a < b]
    for k in range(1, hi + 1):
        cls = [[1 + i * k + c for c in range(k)] for i in range(len(vs))]
        for a, b in es:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return hi + 1


g = degrey.build_G()
pts = list(g.vertices)
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
circle = sorted(g.adj[p_idx])
print(f"pivot {p_idx}, circle of {len(circle)}; "
      f"chi(circle) = {chi(g, circle)}  (bipartite as expected)", flush=True)

# the enriched graph, where every chord point exists, gives the most targets
have, extra = set(pts), []
for i, j in itertools.combinations(circle, 2):
    q = Point(pts[i].x + pts[j].x - pts[p_idx].x,
              pts[i].y + pts[j].y - pts[p_idx].y)
    if q != pts[p_idx] and q not in have:
        have.add(q)
        extra.append(q)
w = build_graph(pts + extra)
pi = w.index_of(pts[p_idx])
circ = w.adj[pi]
aux = [v for v in range(w.n) if v != pi and v not in circ
       and len(w.adj[v] & circ) >= 2]
print(f"enriched: {w.n} vertices, {len(aux)} auxiliaries with two circle "
      f"contacts", flush=True)

# triangles among the auxiliaries: the only shape that could push chi to 5
tri = []
auxset = set(aux)
for a in aux:
    nb = sorted(w.adj[a] & auxset)
    for i in range(len(nb)):
        for j in range(i + 1, len(nb)):
            if nb[j] in w.adj[nb[i]] and a < nb[i]:
                tri.append((a, nb[i], nb[j]))
print(f"  {len(tri)} triangles among them", flush=True)

t0, best = time.time(), 0
for T in tri[:400]:
    c = chi(w, set(circ) | set(T))
    best = max(best, c)
    if c >= 5:
        print(f"  triangle {T}: chi(N(p) u T) = {c}   *** 5-CHROMATIC ON "
              f"{len(circ) + 3} VERTICES ***", flush=True)
        break
print(f"  best chi(N(p) u T) over {min(len(tri), 400)} triangles: {best}  "
      f"[{time.time()-t0:.0f}s]")
print("  (5 would be a 63-vertex 5-chromatic unit-distance graph)")
