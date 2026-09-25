"""Criticality against correlation, measured on both, and they run opposite.

G is nearly 5-vertex-critical: nineteen of thirty probed vertices are
essential, G - v being 4-colourable, one of them verified with an explicit
colouring.  Its correlation at five colours is 1.1.

Sa has 397 vertices where seven suffice for chi = 4, so it should be nowhere
near 4-critical.  Its correlation at four colours is 24.

If that is right then the two quantities run OPPOSITE: the tight graph has no
relation and the loose one has a strong one.  Which would finish the
criticality story off properly, having already withdrawn it once on the
weaker grounds that the pinning graph is the more redundant.

Cheap to settle: Sa is 397 points and a 3-colourability test on it takes
milliseconds.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()


def criticality(name, pts, k, probes=40, seed=555):
    """How often does deleting one vertex drop the chromatic number?"""
    g = build_graph(pts)
    n = g.n
    E = [(min(a, b), max(a, b)) for a, b in g.edges()]
    rng = random.Random(seed)
    sample = rng.sample(range(n), min(probes, n))
    ess = 0
    for v in sample:
        keep = [w for w in range(n) if w != v]
        idx = {w: i for i, w in enumerate(keep)}
        m = len(keep)
        cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
        for a, b in E:
            if a == v or b == v:
                continue
            x, y = idx[a], idx[b]
            for c in range(k):
                cls.append([-(1 + x * k + c), -(1 + y * k + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        if sv.solve():
            ess += 1
        sv.delete()
    print(f"  {name}: {n} points, deleting one vertex makes it "
          f"{k}-colourable in {ess} of {len(sample)} probes  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return ess / len(sample)


# Sa is 4-chromatic: the question is whether Sa - v becomes 3-colourable.
criticality("Sa, dropping to 3", build_Sa(F), 3)
# Y likewise.
criticality("Y, dropping to 3", build_Y(F), 3)
