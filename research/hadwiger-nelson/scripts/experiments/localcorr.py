"""Is correlation local?  If it is, the gadget has to fit in a neighbourhood.

Sa's 576 Moser spindles are seven vertices across and sit inside a disc of
radius about two.  G's single 5-critical subgraph is a thousand vertices and
spans the whole graph.  If the correlation a graph carries is LOCAL -- pairs
at short range, constrained by structure that fits between them -- then that
size difference is the whole story, and no amount of building will help until
a five-colour gadget fits in a neighbourhood.

So measure where the correlated pairs are.  Take Sa's candidates at four
colours -- the pairs differing in every one of twenty-four samples -- and
compare their distance distribution against all pairs of the graph.  If they
sit at short range while the graph's pairs spread out, correlation is local.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def profile(name, pts, k):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    assert sv.solve()
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
    X = np.array([float(p.x) for p in pts])
    Y = np.array([float(p.y) for p in pts])
    cand_d, all_d = [], []
    for i in range(n - 1):
        d = np.sqrt((X[i + 1:] - X[i]) ** 2 + (Y[i + 1:] - Y[i]) ** 2)
        all_d.append(d)
        msk = (C[:, i + 1:] != C[:, i:i + 1]).all(axis=0)
        keep = np.nonzero(msk)[0]
        for off in keep:
            j = int(off) + i + 1
            if (i, j) not in E:
                cand_d.append(float(d[off]))
    allx = np.concatenate(all_d)
    cd = np.array(cand_d)
    if cd.size == 0:
        print(f"  {name}: no candidates", flush=True)
        return
    qs = [50, 75, 90]
    print(f"  {name}: {n} points, {cd.size} correlated pairs of "
          f"{allx.size}; median distance {np.median(cd):.2f} against "
          f"{np.median(allx):.2f} for all pairs; "
          f"fraction within distance 3: {float((cd<3).mean()):.3f} against "
          f"{float((allx<3).mean()):.3f}  [{time.time()-t0:.0f}s]",
          flush=True)
    print(f"    percentiles {qs} of correlated: "
          f"{[round(float(np.percentile(cd,q)),2) for q in qs]}; "
          f"of all pairs: "
          f"{[round(float(np.percentile(allx,q)),2) for q in qs]}",
          flush=True)


profile("Sa at four", build_Sa(F), 4)
profile("Sa at five", build_Sa(F), 5)
profile("G at five", build_G(F, as_graph=False), 5)
