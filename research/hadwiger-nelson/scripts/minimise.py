"""Shrink the 27 classes: which of them are actually needed?

G at five colours stops colouring once twenty-seven of its closable distance
classes are forbidden together, at mean degree 10.9 against its own 10.0 -- a
nine per cent increase in edges, so this is not a density artefact, and it is
the first positive statement about G at five colours in this work: in every
proper 5-colouring, some pair at one of those distances is monochromatic.

Twenty-seven is an artefact of the ORDER though.  They were added smallest
first, which is arbitrary, and a smaller subset may do.  So minimise: drop each
class in turn and keep the drop whenever the rest still fails to colour.  What
survives is irreducible -- every class in it is needed -- and its size is the
honest number.

For scale, the same measure at four colours on Sa is ONE class, D = 16 with
three pairs, which is exactly what de Grey built on.  The gap between the
levels finally has a figure attached, and minimising says how big it really is.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
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
clo = sorted((d for d in byd if closable_distance(d)),
             key=lambda d: len(byd[d]))


def colours(classes):
    cum = list(E)
    for d in classes:
        cum += byd[d]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in cum:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok, len(cum)


cur = clo[:27]
ok, m = colours(cur)
print(f"start: {len(cur)} classes, {m} edges, "
      f"{'colours' if ok else 'does not colour'}  [{time.time()-t0:.0f}s]",
      flush=True)
assert not ok

# drop the LARGEST first: removing a big class is the biggest simplification
order = sorted(cur, key=lambda d: -len(byd[d]))
for d in order:
    trial = [x for x in cur if x != d]
    if not trial:
        break
    ok, m = colours(trial)
    if not ok:
        cur = trial
        print(f"   dropped D={d} ({len(byd[d])} pairs) -> {len(cur)} classes, "
              f"{m} edges, still does not colour  [{time.time()-t0:.0f}s]",
              flush=True)
    else:
        print(f"   D={d} ({len(byd[d])} pairs) is NEEDED"
              f"  [{time.time()-t0:.0f}s]", flush=True)
ok, m = colours(cur)
print(f"\nirreducible: {len(cur)} classes, {m} edges, mean degree "
      f"{2*m/n:.1f}, colours {ok}", flush=True)
print(f"   {[(str(d), len(byd[d])) for d in sorted(cur, key=lambda x: len(byd[x]))]}",
      flush=True)
print(f"   total pairs forbidden: {sum(len(byd[d]) for d in cur)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
