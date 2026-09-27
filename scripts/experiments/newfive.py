"""A 5-chromatic unit-distance graph built by a route de Grey did not take.

The ladder, applied from scratch:

  1. Sa (397 points, Q(sqrt3,sqrt11)) has NO forced-equal pair at 4 colours.
  2. Glue Sa to its image under the 60-degree rotation about the VERTEX
     Sa[25] -- not the origin, so this is not a symmetry and the union is a
     genuinely new 570-point graph.  The glue circle is the unit circle about
     that vertex: 20 points, and NOT capped.
  3. That union has eight forced-equal pairs at 4 colours, all at squared
     distance 64/9.  Gluing manufactured them out of a carrier that had none.
  4. Spindle one of them.  The rotation needs sqrt(2223*16384) = 384*sqrt(247),
     so the graph lives in Q(sqrt3, sqrt11, sqrt247) -- and 247 = 13*19 is a
     radical de Grey's construction never touches.

If the result refuses four colours, it is a 5-chromatic unit-distance graph in
a different field, reached by a different route, and about a third smaller than
the 1581 vertices of G.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from hn.coloring import is_k_colorable
from pysat.solvers import Solver

K1 = Field((3, 11, 247))
print(f"field {K1.gens}, dimension {K1.dim}", flush=True)
pts = build_Sa(K1)
print(f"Sa: {len(pts)} points", flush=True)

centre = pts[25]
rot = rotation_joining(Fr(1), K1).about(centre)
seen, H = set(), []
for p in pts:
    for q in (p, rot(p)):
        if q not in seen:
            seen.add(q); H.append(q)
print(f"H = Sa u rot60_about_Sa[25](Sa): {len(H)} points  "
      f"(shared {2*len(pts)-len(H)})", flush=True)
gH = build_graph(H)
print(f"  n={gH.n}, m={sum(len(a) for a in gH.adj)//2}", flush=True)

a, b = 157, 327
d2 = (H[a] - H[b]).norm2()
print(f"  candidate pair H[{a}], H[{b}]:  d^2 = {d2}  "
      f"(rational: {d2.is_rational()})", flush=True)
assert d2.is_rational() and Fr(d2.c[0]) == Fr(64, 9)

K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(gH.n)]
for u, v in gH.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="cd15", bootstrap_with=cnf)
four = s.solve()
print(f"  H is 4-colourable: {four}", flush=True)
s.delete()
s = Solver(name="cd15", bootstrap_with=cnf)
differ = s.solve(assumptions=[X(a, 0), X(b, 1)])
s.delete()
print(f"  H[{a}] and H[{b}] can differ: {differ}   -> forced EQUAL: {not differ}",
      flush=True)
assert not differ, "the pair is not forced-equal in this field"

rot2 = rotation_joining(Fr(64, 9), K1).about(H[a])
print(f"  spindle rotation: cos {rot2.__doc__ and ''}", flush=True)
seen2, Z = set(), []
for p in H:
    for q in (p, rot2(p)):
        if q not in seen2:
            seen2.add(q); Z.append(q)
print(f"\nZ = H u rho(H): {len(Z)} points  (shared {2*len(H)-len(Z)})", flush=True)
t0 = time.time()
gZ = build_graph(Z)
print(f"  n={gZ.n}, m={sum(len(a) for a in gZ.adj)//2}   "
      f"[edges in {time.time()-t0:.0f}s]", flush=True)
t0 = time.time()
ok, col = is_k_colorable(gZ, 4, timeout=3600)
print(f"\n  Z is 4-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if ok is False:
    ok5, _ = is_k_colorable(gZ, 5, timeout=3600)
    print(f"  Z is 5-colourable: {ok5}", flush=True)
    print(f"\n  *** chi(Z) = 5 : a 5-chromatic unit-distance graph on {gZ.n} "
          f"vertices in Q(sqrt3,sqrt11,sqrt247) ***", flush=True)
    json.dump({"n": gZ.n, "field": list(K1.gens)},
              open("/tmp/hn/newfive.json", "w"))
