"""Does breaking the symmetry raise the correlation?  An anomaly says it might.

Thinning Sa was meant to dilute its Moser spindles and watch the correlation
fall with them.  It mostly did -- 1.45 spindles per point gives ratio 24, 0.98
gives 5.0, 0.70 gives 2.4 -- except at ninety per cent, where the spindle
density FELL to 1.20 and the ratio rose to 251.  Ten times the full graph's.

That is not an artefact: the subgraph is connected, is not 3-colourable, and
returns twenty-four distinct colourings out of twenty-four.  Removing forty
vertices made the colourings ten times more correlated.

A mechanism suggests itself.  Sa is invariant under the twelve-element
dihedral group, so its colourings come in large symmetry orbits and sampling
spreads over them, which decorrelates pairs.  Deleting vertices breaks the
symmetry and the space concentrates.  If that is right it cuts against the
whole symmetrise-then-rotate programme -- G* was built by symmetrising, which
would make things worse, not better.

So: random subgraphs of G at five colours, many of them, and look at the
distribution of the ratio rather than one number.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
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
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    if n < 20:
        return None, n
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return "UNCOLOURABLE", n
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
    return nd / max(ch, 1e-9), n


G = build_G(F, as_graph=False)
base, n0 = ratio(G, 5)
print(f"G whole: {n0} points, ratio {base:.2f}  [{time.time()-t0:.0f}s]",
      flush=True)
rng = random.Random(271828)
best = (0, None)
for trial in range(40):
    frac = rng.choice([0.95, 0.9, 0.85, 0.8, 0.75, 0.7])
    idx = list(range(len(G)))
    rng.shuffle(idx)
    pts = [G[i] for i in sorted(idx[:int(round(frac * len(G)))])]
    r, n = ratio(pts, 5)
    if isinstance(r, str):
        print(f"  trial {trial}: frac {frac}, {n} points -- {r}", flush=True)
        continue
    if r and r > best[0]:
        best = (r, (frac, n))
        print(f"  trial {trial}: frac {frac}, {n} points, ratio {r:.1f}  "
              f"<-- best  [{time.time()-t0:.0f}s]", flush=True)
    elif trial % 10 == 0:
        print(f"  trial {trial}: frac {frac}, {n} points, ratio {r:.1f}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"G at five, 40 random subgraphs: whole {base:.2f}, best "
      f"{best[0]:.1f} at {best[1]}  [{time.time()-t0:.0f}s]", flush=True)
