"""Is the machine list colouring?

Squeezing the pivot's circle into two colours hands every auxiliary point a
LIST: the colours not used by the circle points it is adjacent to. A point
seeing two differently-coloured circle points has a list of k-2; one seeing
two of the same colour has k-1. The minimal refusal of a bipartite
orientation uses both kinds -- six of the first, five of the second, the
latter all at d^2 = 1/3 -- so "confined points against k-2 colours" was never
the right frame. Lists are.

If for every orientation the auxiliary graph fails to be list-colourable with
those lists, the machine is named completely, and an odd cycle is just the
special case where every list is the same pair.
"""
import sys, time, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver
from hn.certify import load_certificate
from hn.graph import build_graph

R = "/home/user/darwin-50/research/hadwiger-nelson"
pts, doc = load_certificate(f"{R}/certificates/pressure3_witness_47.json")
g = build_graph(pts)
piv = 0
circle = set(g.adj[piv])
extras = [v for v in range(g.n) if v != piv and v not in circle]
seen, comps = set(), []
for s in sorted(circle):
    if s in seen:
        continue
    comp, st = [], [s]
    seen.add(s)
    while st:
        x = st.pop()
        comp.append(x)
        for y in g.adj[x] & circle:
            if y not in seen:
                seen.add(y)
                st.append(y)
    comps.append(set(comp))
par, comp_of = {}, {}
for ci, comp in enumerate(comps):
    s = min(comp)
    par[s], front = 0, [s]
    for v in comp:
        comp_of[v] = ci
    while front:
        x = front.pop()
        for y in g.adj[x] & comp:
            if y not in par:
                par[y] = 1 - par[x]
                front.append(y)

sub = build_graph([g.vertices[v] for v in extras])
pos = {v: i for i, v in enumerate(extras)}


def list_colourable(lists):
    var, n = {}, 0
    for i, v in enumerate(extras):
        for c in lists[v]:
            n += 1
            var[(i, c)] = n
    cls = []
    for i, v in enumerate(extras):
        if not lists[v]:
            return False
        cls.append([var[(i, c)] for c in lists[v]])
    for a, b in sub.edges():
        va, vb = extras[a], extras[b]
        for c in set(lists[va]) & set(lists[vb]):
            cls.append([-var[(a, c)], -var[(b, c)]])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


for K in (4, 5):
    t0, bad, sizes = time.time(), 0, collections.Counter()
    for bits in range(1 << len(comps)):
        fixed = {v: par[v] ^ ((bits >> comp_of[v]) & 1) for v in circle}
        lists = {v: [c for c in range(K)
                     if c not in {fixed[u] for u in g.adj[v] & circle}]
                 for v in extras}
        for v in extras:
            sizes[len(lists[v])] += 1
        if not list_colourable(lists):
            bad += 1
    print(f"k={K}: the auxiliary graph is NOT list-colourable in {bad} of "
          f"{1 << len(comps)} orientations; list sizes seen "
          f"{dict(sorted(sizes.items()))}  [{time.time()-t0:.0f}s]", flush=True)
