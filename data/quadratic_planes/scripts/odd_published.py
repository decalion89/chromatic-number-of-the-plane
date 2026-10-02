"""odd_published.py: for each published graph q{d}.json (in the parent directory), whether its direction set U_D
closes a 5-cycle (scan5.has5) and the length of the shortest odd cycle of the graph itself (breadth-first search
from every vertex). On 2 October: the nine graphs with d = 11, 35, 71, 119, 131, 191, 251, 455, 935 have 5-cycles
(among them the three smallest, with 74, 94 and 100 vertices); the other nine have no odd cycle shorter than 7."""
import json, os
from collections import deque
import numpy as np
from units_fast import units_fast
from scan5 import has5

DATA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for f in sorted((f for f in os.listdir(DATA) if f.startswith("q") and f.endswith(".json")), key=lambda s: (len(s), s)):
    g = json.load(open(os.path.join(DATA, f))); d, D = g["d"], g["D"]
    U = np.asarray(units_fast(d, D), dtype=np.int64).reshape(-1, 4)
    n = len(g["points"]); adj = [[] for _ in range(n)]
    for a, b in g["edges"]:
        adj[a].append(b); adj[b].append(a)
    best = None
    for s in range(n):
        dist = [-1] * n; dist[s] = 0; q = deque([s])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if dist[y] < 0:
                    dist[y] = dist[x] + 1; q.append(y)
                elif dist[y] == dist[x]:                 # an edge inside a level closes an odd cycle through s
                    best = 2 * dist[x] + 1 if best is None else min(best, 2 * dist[x] + 1)
    print(f"d={d} D={D}: |U_D| = {len(U)}, 5-cycle among the directions: {has5(U)}; graph: {n} vertices, "
          f"shortest odd cycle {best}", flush=True)
