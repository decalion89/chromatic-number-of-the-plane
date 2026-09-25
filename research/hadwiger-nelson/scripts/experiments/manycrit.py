"""How many 5-critical subgraphs does G actually contain?

The mechanism says correlation needs gadgets at density, and the five-colour
gadget is a 5-critical subgraph.  Sa carries 576 spindles on 397 points -- ten
memberships per point -- and for 500-vertex gadgets to reach the same density
a graph would need about n/50 distinct ones: thirty-two of them inside G.

The question is whether it has one or thirty-two, and it has a cheap side.  If
G - v is 4-COLOURABLE then v lies in every 5-critical subgraph there is, so
there is essentially one; and that case is the fast one, because the solver
only has to exhibit a colouring.  The slow case, G - v still not 4-colourable,
means v is dispensable and there are others.

So budget the query and read the answers that come back quickly.  A budget is
safe here in the direction that matters: it can only fail to find a colouring,
never invent one, so every SAT it reports is a vertex proved essential.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

k, BUDGET = 4, 400000
t0 = time.time()
g = build_G(F)
n = g.n
E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
deg = [len(a) for a in g.adj]
print(f"G: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

rng = random.Random(90210)
probe = rng.sample(range(n), 30)
essential, dispensable, undecided = [], [], []
for nth, v in enumerate(probe):
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
    sv.conf_budget(BUDGET)
    r = sv.solve_limited()
    cf = sv.accum_stats().get("conflicts", 0)
    sv.delete()
    if r is True:
        essential.append(v)
        tag = "4-COLOURABLE: v is in every 5-critical subgraph"
    elif r is False:
        dispensable.append(v)
        tag = "still 5-chromatic without v"
    else:
        undecided.append(v)
        tag = f"undecided within {BUDGET}"
    print(f"  {nth+1}/30  v={v} (degree {deg[v]}): {tag}  ({cf} conflicts)  "
          f"[{time.time()-t0:.0f}s]", flush=True)
print(f"of 30 probed: {len(essential)} essential, {len(dispensable)} "
      f"dispensable, {len(undecided)} undecided  [{time.time()-t0:.0f}s]",
      flush=True)
