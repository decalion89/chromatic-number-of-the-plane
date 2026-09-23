"""Enrich with sigma: every edge grows the two cross-ring points it is missing.

The centre-closure manufactured 1/sqrt3 pairs by the thousand and never helped,
and the reason is now clear: in a triangular lattice all angles are multiples of
30 degrees, while the angle joining the two rings of a hub is
theta0 = arccos(sqrt3/6), which is not one.  So that closure built rings with no
edges to the neighbourhood they were supposed to squeeze.

sigma supplies exactly what was missing.  With sin theta0 = sqrt33/6 the
rotation by theta0 scaled by 1/sqrt3 is exact over Q(sqrt3, sqrt11):

    sigma(x, y) = ( (x - sqrt11 y)/6 , (sqrt11 x + y)/6 ),
    |sigma u| = |u|/sqrt3      and      |u - sigma u| = |u|.

Centre it anywhere: sigma_c(u) = c + sigma(u - c) sits at 1/sqrt3 from c and at
distance |u - c| from u.  So for every EDGE (c,u) of the graph, the two points
c + sigma(u-c) and c + sigmabar(u-c) are new unit-distance neighbours of u lying
on c's 1/sqrt3 ring -- a cross-ring edge, built rather than hoped for.

Iterating that is a genuine enrichment operator, aimed at the one target left:
a vertex that cannot avoid its own 1/sqrt3 ring.  Two things get checked each
round, the second being the whole game.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else 25000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
r11 = F.sqrt(11); sixth = F.rational(Fr(1, 6)); third = F.rational(Fr(1, 3))
def sig(p, c, conj):
    s = -r11 if conj else r11
    dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + (dx - s * dy) * sixth, c.y + (s * dx + dy) * sixth)
pts = {}
for x, y in d["points"]:
    q = Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
    pts[(round(q.fx, 9), round(q.fy, 9))] = q

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
    G = build_graph(list(pts.values())); n = G.n
    E = list(G.edges())
    ok = colourable(G)
    print(f"\nround {rnd}: n={n} edges={len(E)}  5-colourable={ok}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("  *** SIX CHROMATIC ***", flush=True)
        json.dump({"source": NAME, "round": rnd, "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in G.vertices]},
                  open(f"{ROOT}/data/sigma_hit.json", "w"))
        break
    hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
    ring = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
            if abs(dd - 1/3) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == third:
                ring[i].append(j); ring[j].append(i)
    adj = defaultdict(set)
    for x, y in E:
        adj[x].add(y); adj[y].add(x)
    cross = sum(1 for h in ring for v in ring[h] for u in adj[h] if u in adj[v])
    print(f"  rings: {len(ring)} hubs, largest {max((len(v) for v in ring.values()), default=0)};"
          f" cross-ring triangles h-u-v: {cross}   [{time.time()-t0:.0f}s]", flush=True)
    X = lambda v, c: 1 + v * K + c
    base = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in E:
        for c in range(K):
            base.append([-X(x, c), -X(y, c)])
    hit = None
    for h in sorted(ring, key=lambda h: -(len(ring[h]) + len(adj[h])))[:150]:
        cnf = list(base)
        for v in ring[h]:
            for c in range(K):
                cnf.append([-X(h, c), -X(v, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        good = s.solve(); s.delete()
        if not good:
            hit = h; break
    print(f"  ring hub: {hit if hit is not None else 'none of the 150 richest'}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if hit is not None:
        json.dump({"source": NAME, "round": rnd, "hub": hit, "ring": ring[hit]},
                  open(f"{ROOT}/data/sigma_ringhub.json", "w"))
        break
    if rnd == ROUNDS: break
    before = len(pts)
    for c, u in E:
        for (a, b) in ((c, u), (u, c)):
            for conj in (False, True):
                p = sig(G.vertices[b], G.vertices[a], conj)
                k = (round(p.fx, 9), round(p.fy, 9))
                if k not in pts: pts[k] = p
    print(f"  sigma on {len(E)} edges -> +{len(pts)-before} points   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if len(pts) > LIMIT:
        print(f"  stopping at {len(pts)} points", flush=True); break
