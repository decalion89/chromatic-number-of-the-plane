"""What does a K(9/2) map of G actually look like on the plane?

The cyclic-cover reading says such a map is nine distance-avoiding sets
covering every point exactly twice, the pairs being cyclically consecutive.
Nine sets of that kind sitting inside de Grey's graph is a concrete object,
and if they have geometric shape -- stripes, translates, a lattice pattern --
then breaking that shape is a concrete way to attack the rung, rather than
adding points and hoping.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), 9, 2
xy = np.array([[float(pt.x), float(pt.y)] for pt in P])
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
s = Solver(name="cd15", bootstrap_with=cls)
maps = []
for t in range(3):
    if not s.solve():
        break
    mod = s.get_model()
    pos = np.zeros(n, dtype=np.int64)
    for v in range(n):
        for j in range(p):
            if mod[v * p + j] > 0:
                pos[v] = j
                break
    maps.append(pos)
    s.add_clause([-(1 + v * p + int(pos[v])) for v in range(0, n, 37)])
    print(f"map {t}: class sizes "
          f"{[int((pos == j).sum()) for j in range(p)]}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
pos = maps[0]
print("\ngeometry of the nine position classes of map 0:", flush=True)
for j in range(p):
    m = pos == j
    if not m.any():
        print(f"   position {j}: empty")
        continue
    c = xy[m]
    print(f"   position {j}: {m.sum():4d} pts, centroid "
          f"({c[:,0].mean():+.3f},{c[:,1].mean():+.3f}), spread "
          f"{c[:,0].std():.3f}/{c[:,1].std():.3f}", flush=True)
print("\nthe nine distance-avoiding sets (consecutive unions):", flush=True)
adj = collections.defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
for j in range(p):
    J = set(np.nonzero((pos == j) | (pos == (j + 1) % p))[0].tolist())
    bad = sum(1 for v in J for u in adj[v] if u in J) // 2
    print(f"   J_{j} = I_{j} u I_{(j+1)%p}: {len(J):4d} pts, "
          f"{len(J)/n:.4f} of G, {bad} unit pairs inside "
          f"{'(independent as predicted)' if bad == 0 else '(NOT INDEPENDENT)'}",
          flush=True)
if len(maps) > 1:
    agree = [(maps[0] == mm).mean() for mm in maps[1:]]
    print(f"\nlater maps agree with the first on "
          f"{['%.1f%%' % (100*a) for a in agree]} of points", flush=True)
print("DONE", flush=True)
