"""The null control the hill-climb needs.

climbmax reports a union whose best non-edge pair agrees in 25 of 32 sampled
colourings, against 21 of 32 for G alone, and calls that a climb.  Two things
could produce it and they are not the same thing.

 1. The union really constrains the pair, which is the mechanism.
 2. The union is twice the graph.  Twice the vertices is four times the pairs,
    and the maximum of four times as many noisy 32-sample estimates is higher
    for free.  Twice the formula is also a harder search, and a solver under
    5% phase randomisation returns LESS diverse colourings on a bigger
    instance, which raises every pair's agreement at once.

The winning pair sits entirely inside the original copy and the union has ONE
cross edge, so reading (1) into it takes some believing.  The control settles
it: put two copies of G side by side far enough apart that NO cross edge
exists, which is the same n, the same 2m+0 clauses, the same search, and
exactly zero coupling.  Whatever the max reaches there is the free part.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES, k = 32, 5
t0 = time.time()
P = build_G(F, as_graph=False)
g = build_graph(P)
n = g.n
GE = sorted((min(a, b), max(a, b)) for a, b in g.edges())
print(f"G: {n} pts, {len(GE)} edges  [{time.time()-t0:.0f}s]", flush=True)


def sample(m, edges, seed, samples=SAMPLES):
    cls = [[1 + v * k + c for c in range(k)] for v in range(m)]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None
    rng, cols = random.Random(seed), []
    for s in range(samples):
        sv.set_phases([-(1 + w) if rng.random() < .05 else (1 + w)
                       for w in range(m * k)])
        sv.solve()
        mo = sv.get_model()
        cols.append([next(c for c in range(k) if mo[w * k + c] > 0)
                     for w in range(m)])
    sv.delete()
    return np.array(cols, dtype=np.int8)


def topmax(C, edges, m, limit=None):
    """Max agreement over non-edge pairs, restricted to index < limit."""
    es = set(edges)
    hi = limit if limit else m
    top, ti, tj = 0, -1, -1
    for x in range(hi - 1):
        ag = (C[:, x + 1:hi] == C[:, x:x + 1]).sum(axis=0)
        if not len(ag):
            continue
        y = int(ag.argmax())
        if int(ag[y]) > top and (x, y + x + 1) not in es:
            top, ti, tj = int(ag[y]), x, y + x + 1
    return top, ti, tj


# --- G alone, several seeds: how much does the max move on noise? ---
print("\n== G alone, seed to seed ==", flush=True)
for sd in (2718, 31415, 16180, 14142):
    C = sample(n, GE, sd)
    t, i, j = topmax(C, GE, n)
    print(f"  seed {sd}: max {t}/{SAMPLES} at ({i},{j})"
          f"  [{time.time()-t0:.0f}s]", flush=True)

# --- two copies of G, no cross edges at all ---
# Translate by (1000, 0): every cross distance exceeds 900, so no unit edge.
print("\n== G + G translated far away (zero cross edges) ==", flush=True)
DE = GE + [(a + n, b + n) for a, b in GE]
for sd in (2718, 31415):
    C = sample(2 * n, DE, sd)
    t, i, j = topmax(C, DE, 2 * n)
    tin, ii, jj = topmax(C, DE, 2 * n, limit=n)
    print(f"  seed {sd}: max over ALL pairs {t}/{SAMPLES} at ({i},{j});"
          f" max INSIDE copy 1 {tin}/{SAMPLES} at ({ii},{jj})"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("\nDONE", flush=True)
