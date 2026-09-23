"""The free reduction first: a 5-critical graph has minimum degree at least four.

A graph that refuses four colours has a 5-critical subgraph, and in a k-critical
graph every vertex has degree at least k-1.  So any vertex of degree three or
less can be dropped outright -- whatever colouring the rest admits, it has a
free colour left for that vertex -- and the argument iterates, since dropping
one lowers its neighbours' degrees.

That is the 4-core, and it costs no solving at all.  Every solver call in the
shrink hunt is worth two seconds when the vertex is essential and forty-six when
it is removable, so a thousand vertices deleted here is hours not spent there.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
for NAME in sys.argv[1:] or ["five_247.json", "five_247_c.json",
                             "five_dense_2.json", "five_tuned_1_1.json",
                             "five_symmetric.json", "five_twotune_small.json"]:
    try:
        d = json.load(open(f"{ROOT}/data/{NAME}"))
    except FileNotFoundError:
        continue
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    adj = defaultdict(set)
    for x, y in g.edges():
        adj[x].add(y); adj[y].add(x)
    deg = {v: len(adj[v]) for v in range(n)}
    alive = set(range(n))
    q = deque(v for v in alive if deg[v] <= 3)
    while q:
        v = q.popleft()
        if v not in alive: continue
        alive.discard(v)
        for w in adj[v]:
            if w in alive:
                deg[w] -= 1
                if deg[w] == 3: q.append(w)
    m = sum(1 for x, y in g.edges() if x in alive and y in alive)
    print(f"  {NAME:<26s} n={n:<6d} -> 4-core {len(alive):<6d} "
          f"({n-len(alive):>5d} free) edges {m:<6d} degree "
          f"{2*m/max(1,len(alive)):.2f}   [{time.time()-t0:.0f}s]", flush=True)
    if alive and len(alive) < n:
        verts = sorted(alive)
        json.dump({"source": NAME, "n": len(verts),
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in g.vertices[v].x.c],
                               [[t.numerator, t.denominator] for t in g.vertices[v].y.c]]
                              for v in verts]},
                  open(f"{ROOT}/data/core4_{NAME}", "w"))
