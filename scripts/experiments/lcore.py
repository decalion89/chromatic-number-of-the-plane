"""The right object is the whole graph outside the circle, with MIXED lists.

Restricting to the confined set was too narrow.  Squeezing the circle into two
colours gives every other vertex a list whose size depends on what it sees:

    L(v) = 3  if v sees both circle colours   (the confined set)
    L(v) = 4  if v sees exactly one
    L(v) = 5  if v sees none

and the greedy argument completes any list colouring unless some subgraph S
has deg_S(v) >= L(v) for EVERY v in S -- peel vertices with fewer neighbours
than colours and nothing is left.  That maximal surviving subgraph is the
L-CORE, and it is the exact obstruction: empty means the squeeze succeeds and
the pressure is 2, non-empty means it might not.

This is strictly weaker than asking the confined set alone to be 3-degenerate,
because a non-confined neighbour with four colours can still be a neighbour.
So it is the honest version of the condition, and it is what gets measured
here -- on de Grey's Sa, on his Sb and Y, and on the three-hexagon gadget.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from hn.graph import build_graph
from hn.mixed import three_hexagon_gadget


def l_core(g, p, col, k=5):
    """Peel every vertex with fewer neighbours than its list has colours."""
    circle = g.adj[p]
    alive, L = set(), {}
    for v in range(g.n):
        if v == p or v in circle:
            continue
        seen = {col[w] for w in g.adj[v] & circle}
        L[v] = k - len(seen)
        alive.add(v)
    deg = {v: len(g.adj[v] & alive) for v in alive}
    changed = True
    while changed:
        changed = False
        for v in list(alive):
            if deg[v] < L[v]:
                alive.discard(v)
                for u in g.adj[v] & alive:
                    deg[u] -= 1
                changed = True
    return alive, L


def report(name, g, p):
    circle = sorted(g.adj[p])
    cadj = {v: g.adj[v] & set(circle) for v in circle}
    comps, seen = [], set()
    for v in circle:
        if v in seen:
            continue
        stack, part = [v], []
        seen.add(v)
        while stack:
            u = stack.pop()
            part.append(u)
            for w in cadj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(part)
    base = []
    for part in comps:
        col, stack, ok = {}, [(part[0], 0)], True
        while stack:
            u, c = stack.pop()
            if u in col:
                ok &= col[u] == c
                continue
            col[u] = c
            for w in cadj[u]:
                stack.append((w, 1 - c))
        if not ok:
            print(f"{name}: circle component not bipartite -- pressure > 2 free")
            return
        base.append([col, {u: 1 - c for u, c in col.items()}])

    best = 0
    for k in (4, 5):
        worst = None
        for choice in itertools.product(*base):
            col = {}
            for d in choice:
                col.update(d)
            core, L = l_core(g, p, col, k)
            n = len(core)
            if worst is None or n < worst[0]:
                worst = (n, sorted({L[v] for v in core}) if core else [])
        print(f"  {name}, k={k}: circle {len(circle)} in {len(comps)} "
              f"components, {2**len(comps)} orientations; smallest L-core "
              f"{worst[0]}" + (f" (list sizes {worst[1]})" if worst[0] else
                               "  -- EMPTY, so the squeeze always succeeds"),
              flush=True)


for nm, pts in (("de Grey Sa", degrey.build_Sa()),
                ("de Grey Sb", degrey.build_Sb()),
                ("de Grey Y", degrey.build_Y())):
    g = build_graph(pts)
    p = max(range(g.n), key=lambda v: len(g.adj[v]))
    report(nm, g, p)

field, pivot, pts = three_hexagon_gadget()
g = build_graph(pts)
report("three-hexagon gadget", g, g.index_of(pivot))
