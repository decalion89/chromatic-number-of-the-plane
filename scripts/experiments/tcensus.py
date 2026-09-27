"""Rank all 399 translates by their EXACT cross-edge count, and nothing else.

The sampled histogram count is a poor proxy: 451, 445, 439, 435 on the sample
become 1442, 1440, 1442, 1444 exactly, while a lower-ranked one gives 1528 --
and that one costs 17555 conflicts on its dearest pair against the others'
3689.  Five times the cost for six percent more crossing.

So census first, scan second.  Counting only the cross edges, with the copy's
own edges known, is a fraction of a second a translate.


The histogram of t = p - v - q over G puts the best non-trivial translate at
451 cross edges from a sample of 150 of the 1581 points, so the true count
should be near five thousand -- against the 1578 of the best rotation stack
and the six that make de Grey's own step work.

Translates cost nothing arithmetically: no square has to lie in the field, no
rotation has to exist.  They were simply never tried.

Build each union exactly, count what really crosses, and scan every pair at a
closable distance, rerunning anything the budget leaves undecided.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from msqrt import madd, msub, mmul
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

GENS = (3, 5, 7, 11)
DIM = 16
t0 = time.time()
pts = build_G(F, as_graph=False)
P = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c)) for p in pts]
zf = [(float(p.x), float(p.y)) for p in pts]
GE = list(build_graph(pts).edges())
ONE = (Fr(1),) + (Fr(0),) * (DIM - 1)
RAD = []
for m in range(DIM):
    pr = 1
    for b in range(4):
        if m >> b & 1:
            pr *= GENS[b]
    RAD.append(pr ** 0.5)


def fl(c):
    return sum(float(x) * r for x, r in zip(c, RAD))


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


with open("/tmp/hn/gtrans.pkl", "rb") as fh:
    TR = pickle.load(fh)
print(f"G: {len(P)} points; {len(TR)} translates, sampled counts "
      f"{[c for c, _ in TR[:8]]}  [{time.time()-t0:.0f}s]", flush=True)

Pidx = {p: i for i, p in enumerate(P)}
census = []
for ti, (sc, t) in enumerate(TR):
    tz = (fl(t[0]), fl(t[1]))
    img = [(madd(p[0], t[0]), madd(p[1], t[1])) for p in P]
    allp, allz = list(P), list(zf)
    idx = dict(Pidx)
    second = []
    for q, z in zip(img, zf):
        if q in idx:
            second.append(idx[q])
        else:
            idx[q] = len(allp)
            second.append(len(allp))
            allp.append(q)
            allz.append((z[0] + tz[0], z[1] + tz[1]))
    n0, n = len(P), len(allp)
    shared = 2 * n0 - n
    cell = {}
    for i in range(n0):
        a, b = allz[i]
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    E = set()
    for a0, b0 in GE:
        E.add((min(a0, b0), max(a0, b0)))
        x, y = second[a0], second[b0]
        E.add((min(x, y), max(x, y)))
    cross = 0
    for j in range(n0, n):
        a, b = allz[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if abs((a - allz[i][0]) ** 2
                           + (b - allz[i][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = (msub(allp[i][0], allp[j][0]),
                         msub(allp[i][1], allp[j][1]))
                    if norm2(d) == ONE and (i, j) not in E:
                        E.add((i, j))
                        cross += 1
    census.append((cross, shared, n, ti))
    if ti % 40 == 0:
        print(f"  ... {ti}/{len(TR)}, best {max(census)[0]} cross  "
              f"[{time.time()-t0:.0f}s]", flush=True)
census.sort(reverse=True)
print(f"\ntop cross-edge counts: {[c for c, _, _, _ in census[:20]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
print(f"  best: {census[0][0]} cross, {census[0][1]} shared, "
      f"{census[0][2]} points, translate {census[0][3]}", flush=True)
with open("/tmp/hn/tcensus.pkl", "wb") as fh:
    pickle.dump([(c, TR[i][1]) for c, _, _, i in census], fh)
