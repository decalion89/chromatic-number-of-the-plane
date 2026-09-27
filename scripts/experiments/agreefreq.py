"""How OFTEN does the best pair agree?  The graded metric, not the binary one.

Everything here has asked whether a pair agrees in EVERY sampled colouring,
which throws away all the information in between.  The mechanism of Y says
that information is the whole story.  Sa alone does not force (2,0), (-2,0):
its 4-colourings split into those that separate the pair and those that do
not, and the six cross edges of rho kill one class entirely.  If the pair
already agreed in ninety per cent of Sa's colourings, six edges are enough to
finish it; if it agreed at chance, nothing would be.

So measure the FREQUENCY.  For every candidate pair, the fraction of sampled
colourings in which the two ends share a colour, and then the maximum of that
over all pairs -- against 1/k, which is what independence gives.

A base whose best pair sits at 0.9 is one union away from forcing.  One whose
best pair sits at 1/k is not reachable by any union at all.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES = 32
t0 = time.time()


def freq(name, pts, k, seed=2718):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: not colourable", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(seed), []
    for s in range(SAMPLES):
        sv.set_phases([-(1 + w) if rng.random() < .05 else (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    sv.delete()
    C = np.array(cols, dtype=np.int8)
    best, where, hist = 0, None, np.zeros(SAMPLES + 1, dtype=np.int64)
    for i in range(n - 1):
        agree = (C[:, i + 1:] == C[:, i:i + 1]).sum(axis=0)
        np.add.at(hist, agree, 1)
        j = int(agree.argmax())
        if int(agree[j]) > best and (i, j + i + 1) not in E:
            best, where = int(agree[j]), (i, j + i + 1)
    tail = {t: int(hist[t:].sum()) for t in (SAMPLES, SAMPLES - 2,
                                             SAMPLES - 4, SAMPLES // 2)}
    print(f"  {name} at {k}: {n} points; best non-edge pair agrees in "
          f"{best}/{SAMPLES} = {best/SAMPLES:.2f} (chance {1/k:.2f}); pairs "
          f"agreeing in >= t samples: {tail}  [{time.time()-t0:.0f}s]",
          flush=True)


Sa = build_Sa(F)
freq("Sa", Sa, 4)
freq("Sa", Sa, 5)
Y = build_Y(F)
freq("Y", Y, 4)
freq("Y", Y, 5)
G = build_G(F, as_graph=False)
freq("G", G, 5)
rho = rotation_joining(4, F)
Z4 = list({*Sa, *{rho(p) for p in Sa}})
freq("Sa u rho(Sa)", Z4, 4)
