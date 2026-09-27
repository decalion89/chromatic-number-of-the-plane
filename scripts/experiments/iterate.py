"""Iterate the centre operation: G, then G plus every unit triangle's centre.

Adjoining a unit triangle's centre puts a new point at distance exactly 1/sqrt3
from three mutually adjacent old ones -- a whole triangle on its ring at once --
and the centroid of three field points stays in the field, so the iteration
never leaves Q(sqrt3, sqrt11, sqrt247).  Adding points only ever adds
constraints, so chi is non-decreasing along the sequence, and the point set
thickens towards the plane.

One pass on the 803-graph gives 1363 vertices, rings of up to 30 points, and
still a 5-colouring -- and still a 5-colouring even with all 7076 of its 1/sqrt3
pairs forbidden as well.  So iterate, and watch two numbers at each step:

    chi(G_k) <= 5 ?                 the answer itself
    chi(G_k; 1 and 1/sqrt3) <= 5 ?  the distance-forcing statement, which would
                                    say every 5-colouring has a monochromatic
                                    pair at the one radius Niven allows

The second is the weaker and should fall first; it is also the one that, if its
UNSAT core ever concentrates on a single vertex, hands over the hub disjunction
that three rotated copies close.
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
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else 40000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
pts = {}
for x, y in d["points"]:
    q = Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
    pts[(round(q.fx, 9), round(q.fy, 9))] = q
third = F.rational(Fr(1, 3))

def colourable(G, extra=()):
    n = G.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in G.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    for i, j in extra:
        for c in range(K):
            cnf.append([-X(i, c), -X(j, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    r = s.solve(); s.delete()
    return r

for rnd in range(ROUNDS + 1):
    G = build_graph(list(pts.values()))
    n = G.n
    adj = defaultdict(set)
    for x, y in G.edges():
        adj[x].add(y); adj[y].add(x)
    m = sum(len(a) for a in adj.values()) // 2
    ok = colourable(G)
    print(f"\nround {rnd}: n={n} edges={m}  5-colourable={ok}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("  *** SIX CHROMATIC ***", flush=True)
        json.dump({"source": NAME, "round": rnd,
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in G.vertices]},
                  open(f"{ROOT}/data/iterate_hit.json", "w"))
        break
    hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
            if abs(dd - 1/3) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == third:
                pairs.append((i, j))
    r2 = colourable(G, pairs)
    print(f"  1/sqrt3 pairs: {len(pairs)};  chi(G; 1 and 1/sqrt3) <= 5: {r2}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not r2:
        print("  *** THE ONE RADIUS NIVEN ALLOWS IS FORCED ***", flush=True)
        json.dump({"source": NAME, "round": rnd, "pairs": pairs},
                  open(f"{ROOT}/data/forced_root3_radius.json", "w"))
        break
    if rnd == ROUNDS: break
    tri = []
    for u in range(n):
        for v in sorted(adj[u]):
            if v <= u: continue
            for w in sorted(adj[u] & adj[v]):
                if w > v: tri.append((u, v, w))
    before = len(pts)
    for u, v, w in tri:
        c = Point((G.vertices[u].x + G.vertices[v].x + G.vertices[w].x) * third,
                  (G.vertices[u].y + G.vertices[v].y + G.vertices[w].y) * third)
        k = (round(c.fx, 9), round(c.fy, 9))
        if k not in pts: pts[k] = c
    print(f"  {len(tri)} triangles -> +{len(pts)-before} centres   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if len(pts) > LIMIT:
        print(f"  stopping: {len(pts)} points exceeds the limit", flush=True)
        break
