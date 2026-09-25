"""Many hinges, not one.  Does the rigidity come back?

Spindling the lattice once gives a 721-point 4-chromatic graph whose
correlation at four colours is 1.6 -- against Sa's 24, from a base a thousand
times more correlated.  One hinge inherits nothing.

But Sa does not have one hinge.  It carries 1.45 Moser spindles per point,
hundreds of them overlapping, and that density is the one thing separating it
from everything else measured here.  The lattice can be given the same
treatment: its 3-colouring is unique, so every same-class pair is forced, and
every such pair at a closable distance is a hinge waiting to be used.

So take the lattice, find its forced pairs, and union the spindles of many of
them at once.  If the correlation at four colours climbs with the number of
hinges, density is the mechanism and the same thing can be tried at five.  If
it stays at 1.6 however many are added, then a hinge is a hinge and the
rigidity simply does not transfer.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from collections import Counter
from fractions import Fraction as Fr
from itertools import combinations, permutations
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def lattice(side):
    h = F.sqrt(3) / F.rational(2)
    return [Point(F.rational(i) + F.rational(j) / F.rational(2),
                  h * F.rational(j))
            for i in range(-side, side + 1) for j in range(-side, side + 1)]


def measure(name, pts, k):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: NOT {k}-COLOURABLE  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        return None
    rng, cols = random.Random(7919), []
    for s in range(SAMPLES):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    sv.delete()
    C = np.array(cols, dtype=np.int8)
    nd = 0
    for i in range(n - 1):
        nd += int((C[:, i + 1:] != C[:, i:i + 1]).all(axis=0).sum())
    nd -= len(E)
    ch = n * (n - 1) // 2 * ((k - 1) / k) ** SAMPLES
    lut = np.array(list(permutations(range(k))), dtype=np.int8)
    tot, cnt = 0.0, 0
    for a, b in list(combinations(range(SAMPLES), 2))[:30]:
        tot += min(int((C[a] != lut[i][C[b]]).sum())
                   for i in range(len(lut))) / n
        cnt += 1
    r = nd / max(ch, 1e-9)
    print(f"  {name} at {k}: {n} points, {len(E)} edges, ratio {r:.1f}, "
          f"spread {tot/cnt:.3f} (independent {1-1/k:.3f})  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return r


L = lattice(9)
zf = [(float(p.x), float(p.y)) for p in L]
# The unique 3-colouring gives the classes, so the forced pairs are free.
g = build_graph(L)
n0 = g.n
cls = [[1 + v * 3 + c for c in range(3)] for v in range(n0)]
for a, b in g.edges():
    for c in range(3):
        cls.append([-(1 + a * 3 + c), -(1 + b * 3 + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
sv.solve()
m = sv.get_model()
col = [next(c for c in range(3) if m[v * 3 + c] > 0) for v in range(n0)]
sv.delete()

hinges = []
seen = set()
for i in range(n0):
    for j in range(i + 1, n0):
        if col[i] != col[j]:
            continue
        v = (zf[i][0] - zf[j][0]) ** 2 + (zf[i][1] - zf[j][1]) ** 2
        D = Fr(round(v * 55440), 55440)
        if abs(float(D) - v) > 1e-8 or D < Fr(1, 4) or D == 1:
            continue
        if D in seen or not closable_distance(D):
            continue
        seen.add(D)
        hinges.append((D, i, j))
print(f"lattice: {n0} points, {len(hinges)} distinct closable forced "
      f"distances: {[str(d) for d, _, _ in hinges[:10]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

base = set(L)
cur = list(L)
for nth, (D, i, j) in enumerate(hinges[:6]):
    rho = rotation_joining(D, F).about(L[j])
    add = []
    for p in L:
        q = rho(p)
        if q not in base:
            base.add(q)
            add.append(q)
    cur.extend(add)
    measure(f"lattice + {nth+1} hinge(s)", cur, 4)
