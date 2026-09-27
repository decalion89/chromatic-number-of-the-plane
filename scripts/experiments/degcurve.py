"""The degree curve, at four colours and at five, on comparable families.

With the gadget story confounded, what survives is that correlation rises with
the density of constraints -- and two unrelated families already agree on it:
at average degree about 5.6 the triangular lattice gives 1.2 to 2.2 and Sa
thinned to sixty per cent gives 1.7.

So the question worth asking is what the same curve does at five colours.  At
four, degree 9.94 gives 24.0.  At five, G's degree 9.96 gives 1.1.  Thinning G
and measuring across the range settles whether that is a shift of the curve or
a flattening of it: if five colours give about one at EVERY degree, then the
fifth colour does not merely weaken the constraint-density effect, it removes
it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from hn.degrey import build_G, build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def ratio(pts, k):
    g = build_graph(pts)
    n = g.n
    if n < 30:
        return None, 0.0
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return "UNCOLOURABLE", 2 * len(E) / n
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
    return nd / max(ch, 1e-9), 2 * len(E) / n


G = build_G(F, as_graph=False)
Sa = build_Sa(F)
rng = random.Random(104729)
print("G thinned, at five colours:", flush=True)
for frac in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4):
    idx = list(range(len(G)))
    rng.shuffle(idx)
    pts = [G[i] for i in sorted(idx[:int(round(frac * len(G)))])]
    r, deg = ratio(pts, 5)
    rr = r if isinstance(r, str) else f"{r:.1f}"
    print(f"  keep {frac:.1f}: {len(pts)} points, degree {deg:.2f}, "
          f"ratio {rr}  [{time.time()-t0:.0f}s]", flush=True)

print("Sa thinned, at four colours (the comparison):", flush=True)
rng = random.Random(104729)
for frac in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4):
    idx = list(range(len(Sa)))
    rng.shuffle(idx)
    pts = [Sa[i] for i in sorted(idx[:int(round(frac * len(Sa)))])]
    r, deg = ratio(pts, 4)
    rr = r if isinstance(r, str) else f"{r:.1f}"
    print(f"  keep {frac:.1f}: {len(pts)} points, degree {deg:.2f}, "
          f"ratio {rr}  [{time.time()-t0:.0f}s]", flush=True)
