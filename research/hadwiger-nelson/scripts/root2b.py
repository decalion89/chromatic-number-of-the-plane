"""A sqrt2 carrier that actually couples: pick the offset inside the field.

Three translates by w = (1,1) do give every point of G a square centre, but the
copies barely touch -- 2409 vertices carrying only 169 edges between them --
because (1,1) is incommensurate with a triangular lattice: in the basis (1,0),
(1/2, sqrt3/2) its second coordinate is 2/sqrt3.  The offset has to live in the
same field as the carrier, or the translate lands nowhere near it.

Half-offsets with |u|^2 = 1/2 do exist over Q(sqrt3).  Writing u = (p, q) with
p = x + y sqrt3 and q = z + t sqrt3, the conditions are xy + zt = 0 and
x^2 + 3y^2 + z^2 + 3t^2 = 1/2, and (1/4, 1/4, 1/4, -1/4) solves both:

    u = ( (1+sqrt3)/4 , (1-sqrt3)/4 ),   |u|^2 = (4+2sqrt3 + 4-2sqrt3)/16 = 1/2

so 2u has length sqrt2 exactly, inside Q(sqrt3).  Rotating u by multiples of 60
degrees stays in the field and gives more directions to try, and the right one
is whichever puts the most of G a unit away from G + u.

Then G u (G+u) u (G-u) gives h in G the square centre of the diagonal
{h+u, h-u}, with both ends carrying a full copy of the carrier's structure
around them rather than sitting isolated -- which is what a forced colour needs.
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
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
r3 = F.sqrt(3); q4 = F.rational(Fr(1, 4)); half = F.rational(Fr(1, 2))
one = F.rational(1)
u0 = (q4 * (one + r3), q4 * (one - r3))
assert u0[0] * u0[0] + u0[1] * u0[1] == half, "|u|^2 must be 1/2"
cands = [u0]
c60, s60 = half, r3 * half
x, y = u0
for _ in range(5):
    x, y = x * c60 - y * s60, x * s60 + y * c60
    cands.append((x, y))
cands.append((half, half))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
best = None
for k, (ux, uy) in enumerate(cands):
    pts = {}
    for p in P:
        for s in (-1, 0, 1):
            q = Point(p.x + ux * F.rational(s), p.y + uy * F.rational(s))
            pts[(round(q.fx, 9), round(q.fy, 9))] = q
    G = build_graph(list(pts.values()))
    m = sum(1 for _ in G.edges())
    inner = 3 * sum(1 for _ in build_graph(P).edges())
    print(f"  offset {k}: n={G.n} edges={m}  cross edges ~{m-inner}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if best is None or m - inner > best[0]:
        best = (m - inner, k, ux, uy, G, pts)
cross, k, ux, uy, G, pts = best
n = G.n; E = list(G.edges())
print(f"\n  chosen offset {k}: n={n} edges={len(E)}, {cross} cross edges   "
      f"[{time.time()-t0:.0f}s]", flush=True)
idx = {(round(q.fx, 9), round(q.fy, 9)): i for i, q in enumerate(G.vertices)}
trip = []
for p in P:
    h = idx.get((round(p.fx, 9), round(p.fy, 9)))
    a = idx.get((round(float(p.x + ux), 9), round(float(p.y + uy), 9)))
    b = idx.get((round(float(p.x - ux), 9), round(float(p.y - uy), 9)))
    if None not in (a, h, b): trip.append((h, a, b))
adj = defaultdict(set)
for x0, y0 in E:
    adj[x0].add(y0); adj[y0].add(x0)
trip.sort(key=lambda t: -(len(adj[t[0]]) + len(adj[t[1]]) + len(adj[t[2]])))
print(f"  {len(trip)} square-centre triples; richest degrees "
      f"{[len(adj[v]) for v in trip[0]]}   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x0, y0 in E:
    for c in range(K):
        base.append([-X(x0, c), -X(y0, c)])
s = Solver(name="cd19", bootstrap_with=base)
ok = s.solve(); mod = s.get_model() if ok else None
print(f"  5-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if not ok:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"source": NAME, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/root2b_hit.json", "w"))
    sys.exit(0)
col = [next(c for c in range(K) if mod[X(v, c) - 1] > 0) for v in range(n)]
s.delete()
alive = [t for t in trip if col[t[0]] in (col[t[1]], col[t[2]])]
print(f"  one colouring frees {len(trip)-len(alive)} of {len(trip)}; "
      f"{len(alive)} left   [{time.time()-t0:.0f}s]", flush=True)
hit = None
for j, (h, a, b) in enumerate(alive):
    cnf = list(base)
    for c in range(K):
        cnf.append([-X(h, c), -X(a, c)]); cnf.append([-X(h, c), -X(b, c)])
    sv = Solver(name="cd19", bootstrap_with=cnf)
    good = sv.solve(); m2 = sv.get_model() if good else None
    sv.delete()
    if not good:
        hit = (h, a, b); break
    if m2 is not None:
        c2 = [next(c for c in range(K) if m2[X(v, c) - 1] > 0) for v in range(n)]
        alive = [t for t in alive if c2[t[0]] in (c2[t[1]], c2[t[2]])]
    if (j + 1) % 100 == 0:
        print(f"    {j+1} tested, {len(alive)} alive   [{time.time()-t0:.0f}s]",
              flush=True)
print(f"\n  square hub: {hit if hit else 'none'}   [{time.time()-t0:.0f}s]",
      flush=True)
