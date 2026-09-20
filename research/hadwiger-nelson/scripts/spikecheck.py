"""Is the spike the graph, or the sampler?

A ratio of 700 means pairs are surviving twenty-four samples far more often
than chance.  That happens if the graph really constrains them -- or if the
twenty-four colourings are all nearly the same, in which case they carry the
information of one and everything looks correlated.

"Distinct after renaming" does not rule that out: twenty-four colourings
differing in three vertices each are twenty-four distinct tuples.  The test
that does is the distance between them.  Two independent proper k-colourings
of the same graph disagree at about (1 - 1/k) of the vertices; a degenerate
sample disagrees at almost none.

So find a spiking subgraph, and measure how far apart its samples actually
are, against an ordinary one from the same population.
"""
import sys, time, random
from itertools import combinations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES, k = 24, 4
t0 = time.time()
Sa = build_Sa(F)


def sample(pts):
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
        return None, None, None
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
    return nd / max(ch, 1e-9), C, n


def spread(C, n):
    """Mean fraction of vertices on which two samples disagree.

    Colour names are arbitrary, so each pair is compared under the best of the
    k! renamings -- for four colours that is twenty-four permutations, which
    is cheap and removes the only spurious way two identical colourings can
    look far apart.
    """
    from itertools import permutations
    perms = list(permutations(range(k)))
    tot, cnt = 0.0, 0
    for a, b in combinations(range(C.shape[0]), 2):
        best = min(int((C[a] != np.array([p[c] for c in C[b]])).sum())
                   for p in perms)
        tot += best / n
        cnt += 1
    return tot / cnt


rng = random.Random(31337)
found = []
for frac in (0.95, 0.9):
    for trial in range(40):
        idx = list(range(len(Sa)))
        rng.shuffle(idx)
        pts = [Sa[i] for i in sorted(idx[:int(round(frac * len(Sa)))])]
        r, C, n = sample(pts)
        if r is None:
            continue
        if r > 100 or (len(found) < 2 and trial > 35):
            sp = spread(C, n)
            tag = "SPIKE" if r > 100 else "ordinary"
            print(f"  {tag}: frac {frac}, {n} points, ratio {r:.1f}, "
                  f"samples disagree on {sp:.3f} of vertices "
                  f"(independent would be {1-1/k:.3f})  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            found.append((tag, r, sp))
print(f"done: {found}  [{time.time()-t0:.0f}s]", flush=True)
