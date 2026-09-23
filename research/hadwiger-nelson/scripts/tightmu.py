"""mu_5 on the unit circles of the tightest graph: is there a gradient at last?

Every mu_5 measured on a unit circle in this project has read 2 -- in loose
graphs, where any point's neighbourhood can be 2-coloured without the rest of
the graph noticing.  The tight hexagon graph is different in kind: 4159 points
that refuse four and took the solver 452 s to 5-colour, against under a second
for the 803 it was built from.

If its candidate points show mu_5 = 3 anywhere, that is rung one climbed on a
unit circle for the first time -- and it is also a GRADIENT, which is what the
reframing needs.  The target is a graph with essentially one 5-colouring, and
the natural way to walk there is greedy: add the most constrained candidate
point, the one whose neighbourhood is forced to carry the most colours, and
repeat.  Every such point removes colour freedom, and the next candidates' mu
can only rise.  A point with mu_5 = 5 would be a blocked point outright.

So: the richest circle-intersection points of the tight graph, at_most_two on
each, and mu in full wherever that fails.
"""
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
NAME = sys.argv[1] if len(sys.argv) > 1 else "tight_hexagon_4159.json"
TOP = int(sys.argv[2]) if len(sys.argv) > 2 else 150
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
cand = defaultdict(int)
for i in range(n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i: continue
                ex, ey = hx[j] - hx[i], hy[j] - hy[i]
                dd = ex * ex + ey * ey
                if dd <= 1e-12 or dd >= 4.0: continue
                h = math.sqrt(max(0.0, 1.0 - dd / 4.0)); r = math.sqrt(dd)
                mx, my = (hx[i] + hx[j]) / 2, (hy[i] + hy[j]) / 2
                ux, uy = -ey / r, ex / r
                for sg in (1, -1):
                    cand[(round(mx + sg*h*ux, 9), round(my + sg*h*uy, 9))] += 1
existing = {(round(hx[i], 9), round(hy[i], 9)) for i in range(n)}
pool = []
for (kx, ky), c in sorted(cand.items(), key=lambda kv: -kv[1])[:6000]:
    if (kx, ky) in existing: continue
    cx, cy = int(kx // 2), int(ky // 2)
    nb = [i for dx in (-1, 0, 1) for dy in (-1, 0, 1)
          for i in cells.get((cx + dx, cy + dy), ())
          if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
    if len(nb) >= 5: pool.append((len(nb), kx, ky, nb))
pool.sort(key=lambda t: -t[0])
print(f"{NAME}: n={n}; {len(pool)} candidate points with |N| >= 5, richest "
      f"{[t[0] for t in pool[:8]]}   [{time.time()-t0:.0f}s]", flush=True)
ms = MuSolver(G, 5, budget=40_000_000)
print(f"  base colouring {ms.colourable}   [{time.time()-t0:.0f}s]", flush=True)
seen = Counter(); best = None
for m, kx, ky, nb in pool[:TOP]:
    r = ms.at_most_two(nb)
    if r is True:
        seen[2] += 1; continue
    if r is None:
        seen["?"] += 1; continue
    v = ms.mu(nb)
    seen[v] += 1
    print(f"    *** mu_5 = {v} at ({kx:.6f},{ky:.6f}), |N| = {m} ***   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if best is None or (v or 0) > best[0]:
        best = (v, kx, ky, m)
        json.dump({"graph": NAME, "point": [kx, ky], "mu": v, "N": nb},
                  open(f"{ROOT}/data/tight_mu.json", "w"))
print(f"\n  at the {min(TOP,len(pool))} richest points: "
      f"{dict(sorted(seen.items(), key=str))}   [{time.time()-t0:.0f}s]", flush=True)
