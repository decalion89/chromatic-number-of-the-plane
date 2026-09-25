"""Validate the reading of de Grey's last step against his own graph.

The claim: Y carries a pair at distance 4, monochromatic in every proper
4-colouring, and G = Ya u Yb is exactly that pair spindled -- the rotation
closing 4 has sine 3 sqrt7 / 8, which is where sqrt-7 enters his field.

If the claim is right, scanning Y's distance-4 pairs must turn up at least one
that cannot be 2-coloured, and (-2,0) must be one of its ends.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_Y
from hn.geometry import DEGREY_FIELD, Point
from pysat.solvers import Solver

t0 = time.time()
F = DEGREY_FIELD
Y = build_Y(F)
n = len(Y)
zs = [(float(p.x), float(p.y)) for p in Y]
cell = {}
for i, (a, b) in enumerate(zs):
    cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
E = []
for i, (a, b) in enumerate(zs):
    cx, cy = int(a // 1), int(b // 1)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            for j in cell.get((cx + da, cy + db), ()):
                if j <= i:
                    continue
                if abs((a - zs[j][0]) ** 2 + (b - zs[j][1]) ** 2 - 1) > 1e-7:
                    continue
                d = Y[i] - Y[j]
                if d.x * d.x + d.y * d.y == F.rational(1):
                    E.append((i, j))
print(f"Y: {n} vertices, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

P = Point(F.rational(-2), F.zero())
piv = next(i for i, q in enumerate(Y) if q == P)

cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
for a, b in E:
    for c in range(4):
        cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
print(f"  4-colourable? {sv.solve()}  [{time.time()-t0:.0f}s]", flush=True)

four = [(i, j) for i in range(n) for j in range(i + 1, n)
        if abs((zs[i][0] - zs[j][0]) ** 2
               + (zs[i][1] - zs[j][1]) ** 2 - 16) < 1e-7]
print(f"  {len(four)} pairs at distance 4 (float)  [{time.time()-t0:.0f}s]",
      flush=True)
exact = []
for i, j in four:
    d = Y[i] - Y[j]
    if d.x * d.x + d.y * d.y == F.rational(16):
        exact.append((i, j))
print(f"  {len(exact)} confirmed exactly, {sum(1 for i, j in exact if piv in (i, j))}"
      f" of them touching (-2,0)  [{time.time()-t0:.0f}s]", flush=True)

hits = []
for k, (i, j) in enumerate(exact):
    if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
        hits.append((i, j))
        tag = " <- touches (-2,0)" if piv in (i, j) else ""
        print(f"  *** FORCED EQUAL at distance 4: {i},{j}{tag}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    if k and k % 200 == 0:
        print(f"  ... {k}/{len(exact)}  [{time.time()-t0:.0f}s]", flush=True)
print(f"  {len(hits)} forced pairs at distance 4  [{time.time()-t0:.0f}s]",
      flush=True)
if hits:
    print(f"  coordinates of one: {Y[hits[0][0]]}  and  {Y[hits[0][1]]}",
          flush=True)
