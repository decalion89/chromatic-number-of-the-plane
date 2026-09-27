"""Is the tenfold spike a property of ninety per cent, or of one subgraph?

Thinning Sa to ninety per cent of its vertices raised the correlation ratio
from 24 to 251 while lowering the spindle density.  One subgraph is an
anecdote.  If every ninety-per-cent subgraph does it, the effect belongs to
the fraction and there is a mechanism to find; if one in forty does, then some
particular arrangement of 357 points is doing it, and that arrangement is the
only lead in this whole search where correlation went UP unexpectedly.

Forty draws at each of several fractions, same twenty-four samples
throughout, at four colours.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES, k = 24, 4
t0 = time.time()
Sa = build_Sa(F)


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
    return nd / max(ch, 1e-9)


rng = random.Random(31337)
for frac in (0.95, 0.9, 0.85):
    rs = []
    for trial in range(40):
        idx = list(range(len(Sa)))
        rng.shuffle(idx)
        pts = [Sa[i] for i in sorted(idx[:int(round(frac * len(Sa)))])]
        r = ratio(pts)
        if r is not None:
            rs.append(r)
    rs.sort()
    print(f"  frac {frac}: {len(rs)} draws, min {rs[0]:.1f}, median "
          f"{rs[len(rs)//2]:.1f}, max {rs[-1]:.1f}; above 100: "
          f"{sum(1 for x in rs if x > 100)}; above 50: "
          f"{sum(1 for x in rs if x > 50)}  [{time.time()-t0:.0f}s]",
          flush=True)
    print(f"    top five: {[round(x, 1) for x in rs[-5:]]}", flush=True)
