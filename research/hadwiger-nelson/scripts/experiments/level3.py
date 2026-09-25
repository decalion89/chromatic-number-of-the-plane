"""The third level of de Grey's recursion, which needs a field he did not use.

His construction bites twice and the radius doubles: Sa is bitten on its
radius-2 ring by cosine 7/8, and Y is bitten on its radius-4 ring about
(-2, 0) by cosine 31/32, with the second centre a point of the first bitten
ring.  Written as a rule, bite radius rho about c, move to a point of that
ring, bite radius 2*rho about it.

The third step asks for radius 8 and cosine 127/128, whose sine is
sqrt(255)/128 = sqrt(3*5*17)/128.  That number does not exist in
Q(sqrt3, sqrt5, sqrt7, sqrt11), which is why the recursion stops where it
does -- not because the geometry runs out but because the field does.

That is an instruction rather than an obstacle.  Adjoining sqrt17 makes the
turn exact, and the object is small: G has 1581 points and the union has at
most 3162.  Nothing about this is expensive; it just had to be noticed.

The honest expectation is stated in advance.  A bite sharpens a disjunction,
and G has no capped ring at five colours -- the exhaustive scan over all 1581
centres and every radius found none -- so there is nothing here for the third
bite to sharpen.  What the run settles is whether the field was the only thing
in the way.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.field import Field
from hn.degrey import build_G
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
K17 = Field((3, 5, 7, 11, 17))
print(f"field Q(sqrt3, sqrt5, sqrt7, sqrt11, sqrt17), dimension {K17.dim}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
G = build_G(K17, as_graph=False)
print(f"G over the extended field: {len(G)} points  [{time.time()-t0:.0f}s]",
      flush=True)
PIV = Point(K17.rational(-2), K17.zero())
four = K17.rational(16)
ring2 = [p for p in G if (p.x - PIV.x) * (p.x - PIV.x)
         + (p.y - PIV.y) * (p.y - PIV.y) == four]
print(f"the ring of radius 4 about the pivot -- de Grey's second bite -- "
      f"holds {len(ring2)} points  [{time.time()-t0:.0f}s]", flush=True)
if not ring2:
    raise SystemExit("no radius-4 ring: the reading of his construction is wrong")

sixtyfour = K17.rational(64)
c3 = K17.rational(Fr(127, 128))
s3 = K17.sqrt(3) * K17.sqrt(5) * K17.sqrt(17) * K17.rational(Fr(1, 128))
assert c3 * c3 + s3 * s3 == K17.rational(1), "not a rotation"
print(f"the third turn: cos 127/128, sin sqrt(3*5*17)/128 -- exact"
      f"  [{time.time()-t0:.0f}s]", flush=True)

best = None
for i, C in enumerate(ring2):
    on = [p for p in G if (p.x - C.x) * (p.x - C.x)
          + (p.y - C.y) * (p.y - C.y) == sixtyfour]
    if best is None or len(on) > best[1]:
        best = (C, len(on), i)
C, cnt, idx = best
print(f"best third centre: point {idx} of that ring, at "
      f"({float(C.x):.3f},{float(C.y):.3f}); its radius-8 ring holds {cnt} "
      f"points  [{time.time()-t0:.0f}s]", flush=True)
if cnt == 0:
    print("G reaches no point 8 away from any point of the radius-4 ring, so "
          "the third bite has nothing to stitch -- the recursion stops on "
          "geometry as well as on the field", flush=True)
    print("DONE", flush=True)
    raise SystemExit(0)

turn = Rotation(c3, s3).about(C)
seen, P = set(), []
for p in G:
    for q in (p, turn(p)):
        if q not in seen:
            seen.add(q)
            P.append(q)
print(f"union: {len(P)} points  [{time.time()-t0:.0f}s]", flush=True)
g = build_graph(P)
E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
old = set(range(len(G)))
cross = sum(1 for a, b in E if (a in old) != (b in old))
n = g.n
print(f"{n} points, {len(E)} edges, {cross} cross  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + b * k + col)])
for col in range(1, k):
    cls.append([-(1 + col)])
sv = Solver(name="cd15", bootstrap_with=cls)
ok = sv.solve()
sv.delete()
if ok:
    print(f"{k}-colourable  [{time.time()-t0:.0f}s]", flush=True)
else:
    print(f"*** NOT {k}-COLOURABLE -- chi > {k} ***  [{time.time()-t0:.0f}s]",
          flush=True)
    with open("/tmp/claude-0/-home-user-darwin-50/"
              "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
              "WITNESS_level3.pkl", "wb") as f:
        pickle.dump(P, f)
print("DONE", flush=True)
