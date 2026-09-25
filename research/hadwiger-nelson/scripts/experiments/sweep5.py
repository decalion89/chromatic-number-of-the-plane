"""Sweep every distance, not just sqrt3.

The direct analogue of de Grey's lemma fails on G: no vertex has every
5-colouring putting its colour on its sqrt3-sphere.  But sqrt3 was forced at
four colours only because the Moser rotation spindles it, and at five the
right distance need not be the same one.

So sweep.  For each squared distance occurring in G, and each vertex, ask
whether every 5-colouring puts u's colour on that sphere.  The whole sqrt3
pass took eight seconds, so this is affordable.

A hit at distance d would then need a rotation of chord 1/d^2 in the field,
which by the identity 4 - t^2 = (4d^2-1)/d^4 is 3(4d^2-1) being a square
there -- checkable immediately afterwards.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from pysat.solvers import Solver

t0 = time.time()
g = build_G()
pts = list(g.vertices)
E = list(g.edges())
n = len(pts)
print(f"G: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

CUT = 6.0
cells = {}
for i, q in enumerate(pts):
    cells.setdefault((math.floor(q.fx), math.floor(q.fy)), []).append(i)
pairs = Counter()
sphere = {}
rng = range(-int(CUT) - 1, int(CUT) + 2)
for i, q in enumerate(pts):
    cx, cy = math.floor(q.fx), math.floor(q.fy)
    for a in rng:
        for b in rng:
            for j in cells.get((cx + a, cy + b), ()):
                if j <= i:
                    continue
                d2 = (q.fx - pts[j].fx) ** 2 + (q.fy - pts[j].fy) ** 2
                if d2 > CUT * CUT:
                    continue
                k = round(d2 * 1e6)
                pairs[k] += 1
                sphere.setdefault(k, []).append((i, j))
print(f"{len(pairs)} distinct squared distances below {CUT}  "
      f"[{time.time()-t0:.0f}s]", flush=True)


def x(v, c):
    return 1 + v * 5 + c


cls = [[x(v, c) for c in range(5)] for v in range(n)]
for a, b in E:
    for c in range(5):
        cls.append([-x(a, c), -x(b, c)])
solver = Solver(name="cd19", bootstrap_with=cls)

common = [k for k, c in pairs.most_common(60)]
print(f"testing the {len(common)} commonest distances", flush=True)
hits = []
try:
    for k in common:
        star = {}
        for i, j in sphere[k]:
            star.setdefault(i, []).append(j)
            star.setdefault(j, []).append(i)
        best = 0
        for u, S in star.items():
            if len(S) < 4:
                continue
            if not solver.solve(assumptions=[x(u, 0)] + [-x(v, 0) for v in S]):
                hits.append((k / 1e6, u, len(S)))
                print(f"  *** HOLDS at d^2 = {k/1e6:.6f}, vertex {u}, "
                      f"sphere {len(S)}  [{time.time()-t0:.0f}s]", flush=True)
                break
            best = max(best, len(S))
        print(f"  d^2 = {k/1e6:8.4f}: {len(star)} centres, largest sphere "
              f"{best}, no hit  [{time.time()-t0:.0f}s]", flush=True)
finally:
    solver.delete()
print(f"{len(hits)} hits  [{time.time()-t0:.0f}s]", flush=True)
