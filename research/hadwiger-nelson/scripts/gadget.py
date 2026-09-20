"""Does the correlation track the gadget density?  Dilute and see.

The claim is causal: Sa's correlation at four colours comes from its 576 Moser
spindles, and the triangular lattice has none and almost no correlation.  Two
points do not make a curve, and the claim deserves a test that can fail.

Thinning Sa provides one.  A spindle needs all seven of its points, so keeping
a fraction f of the vertices keeps roughly f^7 of the spindles while keeping f
of the points -- the density falls far faster than the graph does.  If the
correlation ratio falls with it, the mechanism is real; if it holds up while
the spindles vanish, the spindles were never what was doing the work.

Measured at a fixed twenty-four samples throughout, at four colours.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES, k = 24, 4
t0 = time.time()
THREE, ONE = F.rational(3), F.rational(1)


def count_spindles(pts):
    n = len(pts)
    zf = [(float(p.x), float(p.y)) for p in pts]
    cell = {}
    for i, (a, b) in enumerate(zf):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    total = 0
    for a in range(n):
        ax, ay = zf[a]
        diag = []
        for cx in range(int((ax - 2.8) // 1), int((ax + 2.8) // 1) + 1):
            for cy in range(int((ay - 2.8) // 1), int((ay + 2.8) // 1) + 1):
                for d in cell.get((cx, cy), ()):
                    if d != a and abs((ax - zf[d][0]) ** 2
                                      + (ay - zf[d][1]) ** 2 - 3.0) < 1e-7:
                        if (pts[a].x - pts[d].x) ** 2 \
                                + (pts[a].y - pts[d].y) ** 2 == THREE:
                            diag.append(d)
        for x in range(len(diag)):
            for y in range(x + 1, len(diag)):
                d, gg = diag[x], diag[y]
                if abs((zf[d][0] - zf[gg][0]) ** 2
                       + (zf[d][1] - zf[gg][1]) ** 2 - 1) > 1e-7:
                    continue
                if (pts[d].x - pts[gg].x) ** 2 \
                        + (pts[d].y - pts[gg].y) ** 2 == ONE:
                    total += 1
    return total


def ratio(pts):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None, n, len(E)
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
    return nd / max(ch, 1e-9), n, len(E)


Sa = build_Sa(F)
rng = random.Random(104729)
for frac in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4):
    idx = list(range(len(Sa)))
    rng.shuffle(idx)
    keep = sorted(idx[:int(round(frac * len(Sa)))])
    pts = [Sa[i] for i in keep]
    sp = count_spindles(pts)
    r, n, m = ratio(pts)
    dens = sp / max(n, 1)
    print(f"  keep {frac:.1f}: {n} points, {m} edges, {sp} spindles "
          f"({dens:.3f}/point), ratio {r if r is None else round(r, 1)}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
