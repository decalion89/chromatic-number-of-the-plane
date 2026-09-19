"""The lists that actually arise, not the ones a chromatic number imagines.

chi(A) = 4 on de Grey's Sa and the pressure at five colours is still 2, so
"the auxiliary graph is 4-chromatic" is not the condition.  The reason is that
the lists are not all equal.  Squeeze the circle into colours {0,1}; then

    an auxiliary seeing BOTH 0 and 1  gets the list {2,3,4}, size 3,
    an auxiliary seeing only one      gets a list of size 4.

Every auxiliary in the first group gets the SAME list, so colouring that group
is ordinary 3-colouring -- and only THAT group has to fail.  The group depends
on which 2-colouring of the circle is chosen, so the condition is:

    pressure(p) > 2 at k = 5  requires the CONFINED set to be 4-chromatic in
    every proper 2-colouring of the circle.

The circle is a union of hexagons, each a 6-cycle with two proper
2-colourings, so the orientations are enumerable exactly.  This measures the
confined set and its chromatic number in each.
"""
import sys, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from hn.graph import build_graph
from pysat.solvers import Solver


def chi(adj, verts, hi=5):
    vs = sorted(verts)
    if not vs:
        return 0
    pos = {v: i for i, v in enumerate(vs)}
    es = [(pos[a], pos[b]) for a in vs for b in adj[a] if b in pos and a < b]
    for k in range(1, hi + 1):
        cls = [[1 + i * k + c for c in range(k)] for i in range(len(vs))]
        for a, b in es:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return hi + 1


g = build_graph(degrey.build_Sa())
p = max(range(g.n), key=lambda v: len(g.adj[v]))
circle = sorted(g.adj[p])
print(f"pivot {p}, circle of {len(circle)}", flush=True)

# proper 2-colourings of the circle's own unit-distance graph
cadj = {v: g.adj[v] & set(circle) for v in circle}
comp, seen = [], set()
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
    comp.append(sorted(part))
print(f"  circle splits into {len(comp)} components of sizes "
      f"{[len(c) for c in comp]}", flush=True)

twocols = []
for c in comp:
    opts, base = [], {}
    stack = [(c[0], 0)]
    ok = True
    while stack:
        u, col = stack.pop()
        if u in base:
            ok &= base[u] == col
            continue
        base[u] = col
        for w in cadj[u]:
            stack.append((w, 1 - col))
    if not ok:
        print(f"  component of size {len(c)} is NOT bipartite: no 2-colouring")
        sys.exit()
    twocols.append([base, {u: 1 - x for u, x in base.items()}])

aux = [v for v in range(g.n)
       if v != p and v not in g.adj[p] and len(g.adj[v] & g.adj[p]) >= 2]
print(f"  {len(aux)} auxiliaries, {2 ** len(comp)} orientations\n", flush=True)

tally = {}
worst = 0
for choice in itertools.product(*twocols):
    col = {}
    for d in choice:
        col.update(d)
    confined = {v for v in aux
                if {col[w] for w in g.adj[v] & g.adj[p]} == {0, 1}}
    c = chi(g.adj, confined)
    tally[c] = tally.get(c, 0) + 1
    worst = max(worst, c)

print("chromatic number of the confined set, over all orientations:")
for c in sorted(tally):
    print(f"   chi = {c}:  {tally[c]} orientations")
print(f"\nmaximum {worst}; pressure > 2 at k = 5 needs 4 in EVERY orientation")
