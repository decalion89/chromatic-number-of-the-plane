"""What the auxiliaries must contain, and what they actually do contain.

Pressure at p is the least number of colours its neighbourhood can be squeezed
into.  Squeezing the circle N(p) into two colours hands every point that sees
two circle neighbours a LIST of the k-2 colours they leave it, and pressure
exceeds 2 exactly when those lists cannot be completed.  So the demand on the
auxiliary graph A -- the points outside the circle with at least two
neighbours on it -- is that it not be (k-2)-choosable:

    k = 4:  lists of size 2.  An odd cycle suffices, which is why de Grey's Sa
            reaches pressure 3 at four colours.
    k = 5:  lists of size 3.  Erdos-Rubin-Taylor: every 2-degenerate graph is
            3-choosable, so A must be at least 3-degenerate, and since it must
            fail 3-choosability it must in particular fail to be 3-COLOURABLE
            somewhere -- A has to contain a 4-chromatic subgraph.

K_4 is not a unit-distance graph, so the cheapest 4-chromatic candidate in the
plane is the Moser spindle.  This measures chi(A) and the degeneracy of A on
every pivot of the structures that reach pressure 3 at k = 4, which is the
sharp test of whether that is the reason they stop at k = 5.
"""
import sys
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.graph import build_graph
from hn.mixed import three_hexagon_gadget
from hn.geometry import Point
from pysat.solvers import Solver


def chi(g_adj, verts, hi=5):
    vs = sorted(verts)
    pos = {v: i for i, v in enumerate(vs)}
    es = [(pos[a], pos[b]) for a in vs for b in g_adj[a] if b in pos and a < b]
    if not vs:
        return 0
    for k in range(1, hi + 1):
        cls = [[1 + i * k + c for c in range(k)] for i in range(len(vs))]
        for a, b in es:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return f">{hi}"


def degeneracy(adj, verts):
    live = {v: len(adj[v] & verts) for v in verts}
    d, rem = 0, dict(live)
    while rem:
        v = min(rem, key=rem.get)
        d = max(d, rem[v])
        for u in adj[v]:
            if u in rem:
                rem[u] -= 1
        del rem[v]
    return d


def report(name, g, pivots):
    print(f"\n{name}: {g.n} vertices", flush=True)
    seen = set()
    for p in pivots:
        circle = g.adj[p]
        aux = {v for v in range(g.n)
               if v != p and v not in circle and len(g.adj[v] & circle) >= 2}
        if not aux or len(aux) in seen:
            continue
        seen.add(len(aux))
        c = chi(g.adj, aux)
        d = degeneracy(g.adj, aux)
        print(f"  pivot {p}: circle {len(circle)}, {len(aux)} auxiliaries, "
              f"chi(A) = {c}, degeneracy {d}"
              + ("   <-- 4-chromatic!" if c not in (0, 1, 2, 3) else ""),
              flush=True)


sa = build_graph(degrey.build_Sa())
deg = sa.degrees()
hubs = sorted(range(sa.n), key=lambda v: -deg[v])[:6]
report("de Grey Sa", sa, hubs)

field, pivot, pts = three_hexagon_gadget()
gad = build_graph(pts)
pi = gad.index_of(pivot)
report("three-hexagon gadget", gad, [pi])

print("\nPREDICTION: chi(A) = 3 everywhere, which is exactly why pressure is 3 "
      "at k = 4 (lists of 2, odd cycle enough) and 2 at k = 5 (lists of 3, "
      "and a 3-colourable A always completes them).")
