"""A carrier with sqrt2 in it, built by three translates.

The unit-square closure needs a vertex h centring a unit square whose diagonal
lies in the graph: h the midpoint of v1, v2 with |v1 - v2| = sqrt2.  No
triangular-lattice carrier has such a pair, since 2 is not Loeschian.  But the
fix is a translation, not a new field.

Take any vector w with |w|^2 = 2 -- Cartesian (1,1) will do, and its
coordinates are rational, so the field never moves -- and form

    G  u  (G + w/2)  u  (G + w).

Then EVERY point p of G supplies the configuration at once: p and p + w are
sqrt2 apart, and their midpoint p + w/2 is a vertex of the middle translate.
Thousands of candidate hubs from one construction, and the translates are
genuinely coupled, since |w/2| = 0.7071 puts many cross pairs at distance
exactly 1.

Two questions, the second being the one that matters:

    is the union still 5-colourable?
    is there a p with c(p + w/2) in {c(p), c(p+w)} forced?

and the second is cheap in the direction it usually answers, since escaping is
satisfiable.  A hit closes the problem with ONE more copy: the 90-degree
rotation about that hub, which is rational.
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
WX = Fr(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
WY = Fr(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
wx, wy = F.rational(WX), F.rational(WY)
half = F.rational(Fr(1, 2))
hx2, hy2 = wx * half, wy * half
two = F.rational(2)
assert (wx * wx + wy * wy) == two, "w must have |w|^2 = 2"
pts = {}
def put(p):
    pts[(round(p.fx, 9), round(p.fy, 9))] = p
for p in P:
    put(p)
    put(Point(p.x + hx2, p.y + hy2))
    put(Point(p.x + wx, p.y + wy))
G = build_graph(list(pts.values())); n = G.n
E = list(G.edges())
print(f"{NAME} + three translates by w=({WX},{WY}), |w|^2=2: n={n} edges={len(E)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
idx = {(round(q.fx, 9), round(q.fy, 9)): i for i, q in enumerate(G.vertices)}
trip = []
for p in P:
    a = idx.get((round(p.fx, 9), round(p.fy, 9)))
    h = idx.get((round(float(p.x + hx2), 9), round(float(p.y + hy2), 9)))
    b = idx.get((round(float(p.x + wx), 9), round(float(p.y + wy), 9)))
    if None not in (a, h, b): trip.append((h, a, b))
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
trip.sort(key=lambda t: -len(adj[t[0]]))
print(f"  {len(trip)} square-centre triples; richest hub degree "
      f"{len(adj[trip[0][0]])}   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])
s = Solver(name="cd19", bootstrap_with=base)
ok = s.solve(); mod = s.get_model() if ok else None
print(f"  5-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if not ok:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"source": NAME, "w": [str(WX), str(WY)],
               "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/root2_hit.json", "w"))
    sys.exit(0)
col = [next(c for c in range(K) if mod[X(v, c) - 1] > 0) for v in range(n)]
s.delete()
alive = [t for t in trip if col[t[0]] in (col[t[1]], col[t[2]])]
print(f"  one exhibited colouring already frees {len(trip)-len(alive)} of "
      f"{len(trip)} triples; {len(alive)} left   [{time.time()-t0:.0f}s]",
      flush=True)
hit = None
for k, (h, a, b) in enumerate(alive):
    cnf = list(base)
    for c in range(K):
        cnf.append([-X(h, c), -X(a, c)])
        cnf.append([-X(h, c), -X(b, c)])
    sv = Solver(name="cd19", bootstrap_with=cnf)
    good = sv.solve()
    m2 = sv.get_model() if good else None
    sv.delete()
    if not good:
        hit = (h, a, b); break
    if m2 is not None:
        c2 = [next(c for c in range(K) if m2[X(v, c) - 1] > 0) for v in range(n)]
        alive = [t for t in alive if c2[t[0]] in (c2[t[1]], c2[t[2]])]
    if (k + 1) % 50 == 0:
        print(f"    {k+1} tested, {len(alive)} still alive   "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\n  square hub: {hit if hit else 'none'}   [{time.time()-t0:.0f}s]",
      flush=True)
if hit:
    json.dump({"source": NAME, "w": [str(WX), str(WY)], "hub": hit[0],
               "diagonal": [hit[1], hit[2]]},
              open(f"{ROOT}/data/square_hub.json", "w"))
