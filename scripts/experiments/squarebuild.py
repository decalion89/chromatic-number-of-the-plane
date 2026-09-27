"""Give the square's centre a neighbourhood: translate a copy onto it.

five_tuned_2_1 carries 17 pairs at sqrt2 -- the distance tuning put them there
on purpose -- and their midpoints are in the field, so adjoining them costs
nothing.  Each midpoint then centres a unit square whose diagonal is already in
the graph, which is exactly what the two-copy closure wants:

    c(h) = c(v1) or c(h) = c(v2) in every 5-colouring
    =>  G u rho_90(G) has no 5-colouring, and rho_90 is rational.

But a bare midpoint has degree one or two, so c(h) is forbidden at most two
colours and escapes trivially.  The closure needs c(N(h)) to cover THREE, since
the two partners can only cover two -- so h needs a real neighbourhood, and the
cheapest real neighbourhood available is the carrier's own.

So translate a copy of the whole graph by t = h - v0, for v0 a vertex of high
degree: h becomes the image of v0 and inherits its 31 neighbours together with
everything that constrains them, while v1 and v2 keep theirs in the original.
The translation is by a field vector, so nothing is extended, and the two copies
overlap wherever they happen to.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_tuned_2_1.json"
HOWMANY = int(sys.argv[2]) if len(sys.argv) > 2 else 3
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
n0 = len(P); two = F.rational(2); half = F.rational(Fr(1, 2))
g0 = build_graph(P)
adj0 = defaultdict(set)
for x, y in g0.edges():
    adj0[x].add(y); adj0[y].add(x)
hx = [q.fx for q in P]; hy = [q.fy for q in P]
cells = defaultdict(list)
for i in range(n0): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
pairs = []
for i in range(n0):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i: continue
                dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
                if abs(dd - 2.0) < 1e-9 and (P[i] - P[j]).norm2() == two:
                    pairs.append((i, j))
pairs.sort(key=lambda p: -(len(adj0[p[0]]) + len(adj0[p[1]])))
v0 = max(range(n0), key=lambda v: len(adj0[v]))
print(f"{NAME} n={n0}: {len(pairs)} sqrt2 pairs; richest diagonal degrees "
      f"{[ (len(adj0[a]), len(adj0[b])) for a,b in pairs[:3] ]}; "
      f"v0 degree {len(adj0[v0])}   [{time.time()-t0:.0f}s]", flush=True)

for idx, (i, j) in enumerate(pairs[:HOWMANY]):
    h = Point((P[i].x + P[j].x) * half, (P[i].y + P[j].y) * half)
    tx, ty = h.x - P[v0].x, h.y - P[v0].y
    pts = {}
    for q in P:
        pts[(round(q.fx, 9), round(q.fy, 9))] = q
        r = Point(q.x + tx, q.y + ty)
        pts[(round(r.fx, 9), round(r.fy, 9))] = r
    G = build_graph(list(pts.values())); n = G.n
    E = list(G.edges())
    adj = defaultdict(set)
    for x, y in E:
        adj[x].add(y); adj[y].add(x)
    key = {(round(q.fx, 9), round(q.fy, 9)): t for t, q in enumerate(G.vertices)}
    H = key[(round(h.fx, 9), round(h.fy, 9))]
    A = key[(round(hx[i], 9), round(hy[i], 9))]
    B = key[(round(hx[j], 9), round(hy[j], 9))]
    print(f"\n  square {idx}: n={n} edges={len(E)}; deg(h)={len(adj[H])}, "
          f"deg(v1)={len(adj[A])}, deg(v2)={len(adj[B])}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    X = lambda v, c: 1 + v * K + c
    base = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in E:
        for c in range(K):
            base.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=base)
    ok = s.solve(); s.delete()
    print(f"    5-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("    *** SIX CHROMATIC ***", flush=True)
        json.dump({"source": NAME, "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in G.vertices]},
                  open(f"{ROOT}/data/squarebuild_hit.json", "w"))
        sys.exit(0)
    cnf = list(base)
    for c in range(K):
        cnf.append([-X(H, c), -X(A, c)]); cnf.append([-X(H, c), -X(B, c)])
    sv = Solver(name="cd19", bootstrap_with=cnf)
    good = sv.solve(); sv.delete()
    print(f"    h can avoid both diagonal ends: {good}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not good:
        print("    *** SQUARE HUB -- two copies close it ***", flush=True)
        json.dump({"source": NAME, "hub": H, "diagonal": [A, B],
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in G.vertices]},
                  open(f"{ROOT}/data/squarebuild_hub.json", "w"))
        break
