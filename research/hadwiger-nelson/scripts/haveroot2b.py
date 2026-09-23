"""Does the graph TUNED to d^2 = 2 actually carry sqrt2 pairs -- and their midpoints?

The claim recorded earlier was that no carrier here has a pair at sqrt2, because
2 is not Loeschian.  That argument is about distances INSIDE one Eisenstein
copy.  Between two copies related by a spindle rotation the squared distance is
a general field element, and this project's distance tuning exists precisely to
choose it: |u + R_phi u|^2 = 2 D^2 (1 + cos phi), with cos phi rational whenever
d^2 and D^2 are.

And d^2 = 2 was one of the four tunings built -- five_tuned_2_1.json, radical
sqrt476 = 4*7*17, 4081 points.  So the question is not whether sqrt2 can be
reached but whether that graph already reaches it, and whether any such pair has
its own MIDPOINT present, which is what the unit-square closure needs:

    h the centre of a unit square whose diagonal {v1,v2} lies in G,
    c(h) = c(v1) or c(h) = c(v2) in every 5-colouring
    =>  G u rho_90(G) has no 5-colouring, and rho_90 is rational.

Asked of every graph in data/, not just the tuned one.
"""
import sys, time, json, os
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
names = sys.argv[1:] or sorted(f for f in os.listdir(f"{ROOT}/data")
                               if f.startswith("five_") and f.endswith(".json"))
for NAME in names:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    if "points" not in d or "field_generators" not in d: continue
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    n = len(P)
    if n > 13000: continue
    hx = [q.fx for q in P]; hy = [q.fy for q in P]
    key = {(round(hx[i], 9), round(hy[i], 9)): i for i in range(n)}
    two = F.rational(2); half = F.rational(Fr(1, 2))
    cells = defaultdict(list)
    for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
    at2 = []
    for i in range(n):
        cx, cy = int(hx[i] // 2), int(hy[i] // 2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cells.get((cx + dx, cy + dy), ()):
                    if j <= i: continue
                    dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
                    if abs(dd - 2.0) < 1e-9 and (P[i] - P[j]).norm2() == two:
                        at2.append((i, j))
    sq = []
    for i, j in at2:
        m = key.get((round((hx[i]+hx[j])/2, 9), round((hy[i]+hy[j])/2, 9)))
        if m is not None: sq.append((m, i, j))
    tag = "  <-- SQUARE CENTRES" if sq else ""
    print(f"  {NAME:<28s} n={n:<6d} pairs at sqrt2: {len(at2):<6d} "
          f"square centres: {len(sq)}{tag}   [{time.time()-t0:.0f}s]", flush=True)
    if sq:
        g = build_graph(P)
        adj = defaultdict(set)
        for x, y in g.edges(): adj[x].add(y); adj[y].add(x)
        sq.sort(key=lambda t: -(len(adj[t[0]]) + len(adj[t[1]]) + len(adj[t[2]])))
        print(f"    richest: hub deg {len(adj[sq[0][0]])}, diagonal "
              f"{len(adj[sq[0][1]])},{len(adj[sq[0][2]])}", flush=True)
        json.dump({"graph": NAME, "squares": [[a, b, c] for a, b, c in sq[:2000]]},
                  open(f"{ROOT}/data/squares_{NAME}", "w"))
