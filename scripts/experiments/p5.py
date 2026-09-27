"""How near does G come to closing a directed path on five vertices?

chi_c(G) >= 5 iff every acyclic orientation has a directed 4-edge path whose
ends are adjacent.  G is 5-chromatic so Gallai-Roy already gives every
orientation a directed 4-edge path; the ends are just never a unit apart.
So the question is quantitative: in G's best orientation, how far apart ARE
the ends of those paths, and does the distribution come near 1?

A spike just off 1 would say the graph is close and a small addition could
close it.  A distribution with a hole at 1 would say the geometry keeps the
ends away on purpose, which is a different and more interesting obstruction.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, collections
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), 9, 2
xy = np.array([[float(t.x), float(t.y)] for t in P])
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve()
mod = s.get_model()
pos = np.zeros(n, dtype=np.int64)
for v in range(n):
    for j in range(p):
        if mod[v * p + j] > 0:
            pos[v] = j
            break
out = collections.defaultdict(list)
for a, c in E:
    if pos[a] < pos[c]:
        out[a].append(c)
    else:
        out[c].append(a)
print(f"G: {n} pts, {len(E)} edges; orientation by position, "
      f"{sum(len(v) for v in out.values())} arcs  [{time.time()-t0:.0f}s]",
      flush=True)
random.seed(4)
ends, seen = [], set()
starts = list(out)
random.shuffle(starts)
for a in starts[:600]:
    for bb in out[a]:
        for c in out.get(bb, []):
            for d in out.get(c, []):
                for e in out.get(d, []):
                    if (a, e) in seen:
                        continue
                    seen.add((a, e))
                    ends.append(float(np.linalg.norm(xy[a] - xy[e])))
ends = np.array(ends)
print(f"{len(ends)} distinct end pairs of directed 4-edge paths"
      f"  [{time.time()-t0:.0f}s]", flush=True)
if len(ends):
    print(f"   distance range {ends.min():.4f} .. {ends.max():.4f}, "
          f"median {np.median(ends):.4f}", flush=True)
    for lo, hi in ((0.0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.01),
                   (1.01, 1.1), (1.1, 1.5), (1.5, 2.0), (2.0, 9.0)):
        c0 = int(((ends >= lo) & (ends < hi)).sum())
        print(f"   [{lo:.2f},{hi:.2f}): {c0:6d}  {100*c0/len(ends):5.2f}%",
              flush=True)
    near = ends[(ends > 0.9) & (ends < 1.1)]
    print(f"\n   {len(near)} pairs within 10% of a unit apart; closest to 1: "
          f"{sorted(abs(near - 1))[:6] if len(near) else 'none'}", flush=True)
print("DONE", flush=True)
