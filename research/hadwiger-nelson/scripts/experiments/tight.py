"""Which EXISTING vertices have tight neighbourhoods, and what do they look like?

The bipartite theorem closes the single-point attack and points at its own
corollary.  If H is 6-chromatic and vertex-critical then for every vertex v,
H - v is 5-colourable and N(v) uses all five colours in every 5-colouring of
it -- and N(v) is bipartite.  So a 6-chromatic unit-distance graph must carry,
at every vertex, a BIPARTITE neighbourhood that the rest of the graph forces
to five colours.

New points never got past two.  The question is whether existing vertices do
better, because the same measurement applies to them: for a vertex v of G, the
minimum over proper 5-colourings of the number of colours on N(v).  A vertex
whose minimum is 4 is as tight as anything can be at five colours -- its
neighbourhood is one colour from rainbow -- and its structure is the template
for what a forced neighbourhood looks like.

Reported per vertex: degree, and the minimum colours its neighbourhood is
forced to.  If the distribution reaches 4 the template exists and is worth
copying; if it also stops at 2, the flexibility is the graph's and not the
candidates'.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from collections import Counter, defaultdict
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()


def run(name, P, kk=k):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    n = len(P)
    adj = defaultdict(list)
    for a, b in E:
        adj[a].append(b)
        adj[b].append(a)
    cls = [[1 + v * kk + c for c in range(kk)] for v in range(n)]
    for a, b in E:
        for c in range(kk):
            cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    assert sv.solve(), f"{name} is not {kk}-colourable"
    dist = Counter()
    bydeg = defaultdict(Counter)
    best = (0, None)
    for v in range(n):
        nb = adj[v]
        if len(nb) < 2:
            continue
        lo = kk
        for j in range(1, kk + 1):
            ass = [-(1 + u * kk + c) for u in nb for c in range(j, kk)]
            if sv.solve(assumptions=ass):
                lo = j
                break
        dist[lo] += 1
        bydeg[len(nb)][lo] += 1
        if lo > best[0]:
            best = (lo, v, len(nb))
    sv.delete()
    print(f"\n{name} at {kk} colours: {n} pts, {len(E)} edges", flush=True)
    print(f"   minimum colours forced on N(v), over the vertices: "
          f"{dict(sorted(dist.items()))}  [{time.time()-t0:.0f}s]", flush=True)
    tops = sorted(bydeg)[-5:]
    for d in tops:
        print(f"      degree {d:3d}: {dict(sorted(bydeg[d].items()))}",
              flush=True)
    if best[0] >= 2:
        print(f"   tightest: vertex {best[1]} of degree {best[2]}, "
              f"neighbourhood forced to {best[0]} colours", flush=True)


run("G", build_G(K, as_graph=False), 5)
run("Sa", build_Sa(K), 4)
run("Y", build_Y(K), 4)
print("\nDONE", flush=True)
