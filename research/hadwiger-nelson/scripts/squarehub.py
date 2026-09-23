"""Adjoin the midpoints: the tuned-to-two graph already carries the diagonals.

The recorded claim was that no carrier here has a pair at sqrt2, and it was too
strong.  The Loeschian argument -- a^2 + ab + b^2 never equals 2 -- governs
distances INSIDE one Eisenstein copy.  The distance tuning crosses copies, and
|u + R_phi u|^2 = 2 D^2 (1 + cos phi) is exactly the knob that chooses the
answer; d^2 = 2 was one of the four tunings built.  Asked of every graph in
data/: five_tuned_2_1 has 17 pairs at sqrt2, five_tuned_4_1 has 8, and every
other graph has none, which is why the check on the sigma-enriched 803 came
back zero and got over-generalised.

Their midpoints are not vertices, but a midpoint is (u+v)/2 and the field is
closed under that, so adjoining them costs nothing and no extension.  Each one
then centres a unit square whose diagonal is already in the graph -- and the
closure needs only

    c(h) = c(v1)  or  c(h) = c(v2)   in every 5-colouring,

after which G union its 90-degree rotation about h has no 5-colouring.  Two
copies, and the rotation is rational.

Three things to read: whether the enriched graph still admits five colours, how
rich each new hub's own neighbourhood is -- c(h) is only cornered if c(N(h))
already covers three colours -- and whether any hub is cornered.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_tuned_2_1.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
n0 = len(P); two = F.rational(2); half = F.rational(Fr(1, 2))
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
print(f"{NAME} n={n0}: {len(pairs)} pairs at sqrt2   [{time.time()-t0:.0f}s]",
      flush=True)
pts = {(round(q.fx, 9), round(q.fy, 9)): q for q in P}
mids = []
for i, j in pairs:
    m = Point((P[i].x + P[j].x) * half, (P[i].y + P[j].y) * half)
    k = (round(m.fx, 9), round(m.fy, 9))
    if k not in pts: pts[k] = m
    mids.append((k, i, j))
G = build_graph(list(pts.values())); n = G.n
E = list(G.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
key = {(round(q.fx, 9), round(q.fy, 9)): t for t, q in enumerate(G.vertices)}
trips = []
for k, i, j in mids:
    h = key[k]
    a = key[(round(hx[i], 9), round(hy[i], 9))]
    b = key[(round(hx[j], 9), round(hy[j], 9))]
    trips.append((h, a, b))
print(f"  + {n - n0} midpoints -> n={n} edges={len(E)}   [{time.time()-t0:.0f}s]",
      flush=True)
print(f"  hub degrees: {sorted((len(adj[h]) for h, _, _ in trips), reverse=True)}",
      flush=True)
print(f"  diagonal degrees: "
      f"{sorted(((len(adj[a]), len(adj[b])) for _, a, b in trips), reverse=True)[:4]}",
      flush=True)
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
    json.dump({"source": NAME, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/square_hit.json", "w"))
    sys.exit(0)
col = [next(c for c in range(K) if mod[X(v, c) - 1] > 0) for v in range(n)]
s.delete()
alive = [t for t in trips if col[t[0]] in (col[t[1]], col[t[2]])]
print(f"  one colouring frees {len(trips)-len(alive)} of {len(trips)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
for h, a, b in alive:
    cnf = list(base)
    for c in range(K):
        cnf.append([-X(h, c), -X(a, c)]); cnf.append([-X(h, c), -X(b, c)])
    sv = Solver(name="cd19", bootstrap_with=cnf)
    good = sv.solve(); sv.delete()
    if not good:
        print(f"  *** SQUARE HUB h={h}, diagonal ({a},{b}) ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        json.dump({"source": NAME, "hub": h, "diagonal": [a, b]},
                  open(f"{ROOT}/data/square_hub.json", "w"))
        break
else:
    print(f"  no square hub among {len(alive)}   [{time.time()-t0:.0f}s]",
          flush=True)
