"""The one point that can be put in the gap: G* at five colours.

At four colours, Sa has 397 vertices where the minimum for chi = 4 is seven
-- fifty-seven times over -- and its colour relation sits twenty times above
chance.  At five colours, G has 1581 where about five hundred suffice, three
times over, and its relation sits AT chance and verifies to nothing.

Between three times the minimum and fifty-seven times, nothing has been
measured at five colours.  G* is the only object here that falls inside it:
13873 vertices against the same five hundred, twenty-eight times over.  If
size relative to the minimum is what turns the relation on, G* is where it
should start to show; if G* also sits at chance, then size is not it either
and five colours is simply out of reach of this whole family.

Sixty samples rather than twenty-four, because ninety-six million pairs need
(4/5)^60 to cut them to something a solver can check one by one.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.graph import build_graph
from hn.homcol import agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/hn/")
SAMPLES, k = 60, 5
t0 = time.time()
PIV = Point(F.rational(-2), F.zero())
rot = _rot60(F).about(PIV)
G = build_G(F, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
g = build_graph(Gs)
n = g.n
E = set((min(a, b), max(a, b)) for a, b in g.edges())
print(f"G*: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in sorted(E):
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="g4", bootstrap_with=cls)
assert sv.solve()
rng, cols = random.Random(2027), []
for s in range(SAMPLES):
    sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                   for w in range(n * k)])
    sv.solve()
    m = sv.get_model()
    cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                 for w in range(n)])
    if (s + 1) % 10 == 0:
        print(f"  ... {s+1}/{SAMPLES} samples  [{time.time()-t0:.0f}s]",
              flush=True)

same = [(i, j) for i, j in agreeing_pairs(cols, colours=k)
        if (i, j) not in E]
C = np.array(cols, dtype=np.int8)
diff = []
for i in range(n - 1):
    m = (C[:, i + 1:] != C[:, i:i + 1]).all(axis=0)
    for off in np.nonzero(m)[0]:
        j = int(off) + i + 1
        if (i, j) not in E:
            diff.append((i, j))
pairs = n * (n - 1) // 2
chance = pairs * ((k - 1) / k) ** SAMPLES
print(f"G* at five: {pairs} pairs, chance {chance:.1f}, "
      f"{len(same)} forced-same candidates, {len(diff)} forced-different "
      f"candidates, ratio {len(diff)/max(chance,1e-9):.2f}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

fs = [(i, j) for i, j in same
      if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
fd = [(i, j) for i, j in diff[:4000]
      if not sv.solve(assumptions=[1 + i * k, 1 + j * k])]
sv.delete()
print(f"VERIFIED: {len(fs)} forced-same, {len(fd)} forced-different beyond "
      f"the edges (of {min(len(diff), 4000)} checked)  "
      f"[{time.time()-t0:.0f}s]", flush=True)
if fs or fd:
    with open(SC + "relbig.pkl", "wb") as fh:
        pickle.dump(([(str(p.x), str(p.y)) for p in Gs], fs, fd), fh)
    print("  *** G* HAS A NON-TRIVIAL COLOUR RELATION AT FIVE COLOURS ***",
          flush=True)
