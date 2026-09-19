"""The lift, asked as the right question: choosability, over every orientation.

The circle uses two colours, so an auxiliary seeing both is left with exactly
the other k-2 -- and every such point gets the SAME list. At five colours that
makes the confined subgraph's requirement plain: it must be 4-chromatic. The
points seeing only one circle colour get k-1 and can only help by propagation.

So: take the pivot's circle, generate every exact point one away from two of
its points, and for each orientation ask both questions -- is the confined
subgraph 4-chromatic, and does the whole list instance refuse. The worst
orientation is the answer.
"""
import sys, time, collections, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.forced import circle_hexagons, list_colourable
from hn.graph import build_graph
from hn.mixed import circle_intersections

R = "/home/user/darwin-50/research/hadwiger-nelson"
pts, doc = load_certificate(f"{R}/certificates/pressure3_witness_47.json")
g0 = build_graph(pts)
field = pts[0].x.field
one = field.rational(1)
piv0 = 0
circle0 = sorted(g0.adj[piv0])

t0 = time.time()
cand, seenq = [], set(pts)
for a, b in itertools.combinations(circle0, 2):
    if float(g0.vertices[a].dist2(g0.vertices[b])) > 3.999:
        continue
    for q in circle_intersections(g0.vertices[a], one, g0.vertices[b], one):
        if q not in seenq:
            seenq.add(q)
            cand.append(q)
print(f"{len(cand)} candidate auxiliaries  [{time.time()-t0:.0f}s]", flush=True)

allpts = list(pts) + cand
G = build_graph(allpts)
piv = G.vertices.index(g0.vertices[piv0])
circle = set(G.adj[piv])
aux = [v for v in range(G.n) if v != piv and v not in circle]
comps, par, comp_of = circle_hexagons(G, piv)
print(f"combined {G}   circle {len(circle)} in {len(comps)} components of "
      f"sizes {sorted(len(c) for c in comps)}, {len(aux)} auxiliaries",
      flush=True)

K = 5
worst = None
for bits in range(1 << len(comps)):
    fixed = {v: par[v] ^ ((bits >> comp_of[v]) & 1) for v in circle}
    lists, conf = {}, []
    for v in aux:
        used = {fixed[u] for u in G.adj[v] & circle}
        lists[v] = [c for c in range(K) if c not in used]
        if len(used) == 2:
            conf.append(v)
    sub = build_graph([G.vertices[v] for v in conf])
    chi = 2
    while chi < 6 and not is_k_colorable(sub, chi)[0]:
        chi += 1
    refused = not list_colourable(G, aux, lists)
    if worst is None or (chi, refused) < worst[:2]:
        worst = (chi, refused, bits, len(conf))
        print(f"  orientation {bits:0{len(comps)}b}: {len(conf)} confined, "
              f"their chi = {chi}, whole instance refused: {refused}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"WORST: confined chi {worst[0]}, refused {worst[1]}, "
      f"{worst[3]} confined.  k=5 needs chi >= 4 or a refusal everywhere  "
      f"[{time.time()-t0:.0f}s]", flush=True)
