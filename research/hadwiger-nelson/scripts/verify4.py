"""Verify the four-colour claim that the whole contrast rests on.

Sa at four colours leaves 1548 pairs differing in every one of twenty-four
samples, against 79 expected by chance -- nineteen times over -- and that
ratio is the anchor the five-colour negatives are read against.  But not one
of those 1548 was ever put to the solver: at four colours a forced-different
proof runs to minutes, so the scan reported counts and moved on.

Counts are not a claim.  If none of them is genuinely forced, the ratio is
sampling bias and the contrast evaporates.  So take a sample of them, check
each with no budget anywhere, and see.

Small graphs go first as a sanity check, where the answer is known: a rhombus
at three colours forces its apexes to agree -- a constrained pair that is not
an edge -- while K4 at four colours has no non-edges at all, so criticality by
itself guarantees nothing.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()


def relation_of(name, n, E, k, pairs=None, cap=None):
    Es = set((min(a, b), max(a, b)) for a, b in E)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(Es):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    assert sv.solve(), f"{name} is not {k}-colourable"
    if pairs is None:
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)
                 if (i, j) not in Es]
    if cap:
        pairs = pairs[:cap]
    same = [(i, j) for i, j in pairs
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    diff = [(i, j) for i, j in pairs
            if not sv.solve(assumptions=[1 + i * k, 1 + j * k])]
    sv.delete()
    print(f"  {name}: {len(pairs)} non-edge pairs checked -> {len(same)} "
          f"forced-same, {len(diff)} forced-different  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return same, diff


# A rhombus at three colours: the apexes are a constrained non-edge.
relation_of("rhombus at 3", 4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)], 3)
# K4 at four colours: 4-critical, and no non-edges at all.
relation_of("K4 at 4", 4, [(a, b) for a in range(4) for b in range(a + 1, 4)],
            4)
# The Moser spindle: 7 vertices, 4-critical, eleven edges.
MOSER = [(0, 1), (0, 2), (1, 2), (0, 3), (0, 4), (3, 4), (2, 5), (1, 6),
         (5, 6), (3, 5), (4, 6)]
relation_of("Moser spindle at 4", 7, MOSER, 4)

# And the one that matters: Sa's candidates at four colours, sampled.
Sa = build_Sa(F)
g = build_graph(Sa)
n, k = g.n, 4
E = set((min(a, b), max(a, b)) for a, b in g.edges())
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, b in sorted(E):
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="g4", bootstrap_with=cls)
sv.solve()
rng, cols = random.Random(6180), []
for s in range(24):
    sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                   for w in range(n * k)])
    sv.solve()
    m = sv.get_model()
    cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                 for w in range(n)])
sv.delete()
C = np.array(cols, dtype=np.int8)
cand = []
for i in range(n - 1):
    msk = (C[:, i + 1:] != C[:, i:i + 1]).all(axis=0)
    for off in np.nonzero(msk)[0]:
        j = int(off) + i + 1
        if (i, j) not in E:
            cand.append((i, j))
print(f"Sa at four: {len(cand)} candidates from 24 samples "
      f"(chance {n*(n-1)//2*(3/4)**24:.0f})  [{time.time()-t0:.0f}s]",
      flush=True)
rng.shuffle(cand)
relation_of("Sa at 4, 40 sampled candidates", n, sorted(E), 4,
            pairs=cand[:40])
