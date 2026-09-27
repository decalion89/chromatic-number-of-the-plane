"""Minimise the obstruction pair by pair, not class by class.

The five classes are irreducible AS CLASSES -- drop any one and G colours
again -- but that says nothing about the pairs inside them.  144 pairs at
distance 7 were kept or dropped together, and there is no reason the whole 144
are needed.

So minimise properly: run through the 229 pairs, drop each one, and keep the
drop whenever the formula still fails to colour.  What survives is irreducible
pair by pair -- every single pair in it is load-bearing.

The shape of the survivor is what matters.  If it collapses onto ONE distance
class, the statement becomes "in every 5-colouring some pair at distance d is
monochromatic", which is exactly de Grey's hypothesis at five colours and the
ordinary two-copy spindle consumes it.  If it stays spread across five
distances, a multispindle is unavoidable and this says how many centres and
angles it needs.

Largest class first, since that is where the slack most likely is.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
KEEP = [Fr(7), Fr(9), Fr(17, 2), Fr(16), Fr(15, 16)]
P = build_G(K, as_graph=False)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
n = len(P)
E = sorted(set((min(a, b), max(a, b))
               for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
byd = defaultdict(list)
for i in range(n - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in Eset:
            v = int(sq[off, 0])
            if v:
                byd[Fr(v, D2)].append((i, j))

cur = [(D, p) for D in KEEP for p in byd[D]]
print(f"start: {len(cur)} pairs over {len(KEEP)} classes"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def colours(pairs):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    for _, (a, b) in pairs:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


assert not colours(cur)
# Iterate over a fixed SNAPSHOT, by identity.  The first version computed the
# order from the original list and then indexed into the SHRINKING one, so
# after each drop the indices pointed at different elements and some pairs
# were never tested at all -- which makes "irreducible" unearned.  Passes
# repeat until a whole sweep drops nothing, since removing one pair can free
# another.
sweep = 0
while True:
    sweep += 1
    snapshot = sorted(cur, key=lambda x: -float(x[0]))
    dropped = 0
    for target in snapshot:
        if target not in cur:
            continue
        trial = [x for x in cur if x is not target]
        if colours(trial):
            continue
        cur = trial
        dropped += 1
        if dropped % 25 == 0:
            c = Counter(str(D) for D, _ in cur)
            print(f"   sweep {sweep}, dropped {dropped}, {len(cur)} left: "
                  f"{dict(c)}  [{time.time()-t0:.0f}s]", flush=True)
    c = Counter(str(D) for D, _ in cur)
    print(f"  sweep {sweep}: dropped {dropped}, {len(cur)} left, "
          f"{dict(sorted(c.items()))}  [{time.time()-t0:.0f}s]", flush=True)
    if dropped == 0:
        break

c = Counter(str(D) for D, _ in cur)
print(f"\nirreducible pair by pair after {sweep} sweeps: "
      f"{len(cur)} pairs over {len(c)} classes",
      flush=True)
print(f"   {dict(sorted(c.items()))}  [{time.time()-t0:.0f}s]", flush=True)
if len(c) == 1:
    print("   >>> ONE distance class -- this is de Grey's hypothesis at five",
          flush=True)
for D, (a, b) in sorted(cur, key=lambda x: str(x[0]))[:12]:
    print(f"      D={D}: {a},{b}", flush=True)
print("DONE", flush=True)
