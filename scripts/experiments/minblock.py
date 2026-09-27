"""How few distances can carry the disjunction?

S is a set in G whose palette is 4 of 5, so every 5-colouring gives it a
monochromatic class of at least ceil(|S|/4) points, and that class is
independent.  Forbidding every rational distance at once leaves no
independent class of that size, so in every colouring some monochromatic pair
inside S sits at one of those distances -- a disjunction over a named, finite
list rather than over nothing.

The list is what matters now.  A disjunction over many distances is hard to
use; over one distance it is a forced pair, and a forced pair at a closable
distance spindles into a sixth colour.  So: shrink the list.  Drop a distance,
ask whether an independent class of the forced size can now escape, and keep
the drop only when it cannot.  What survives is a minimal blocking set -- no
proper subset of it carries the disjunction.

Every step is a finite satisfiability question over |S| variables, so the
answer is exact.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound, closable_distance
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t0 = time.time()
SC = ("/tmp/hn/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
Eset = set(E)
n = len(P)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
S, _ = pickle.load(open(SC + f"grown_{k}_{TARGET}.pkl", "rb"))
S = sorted(S)
m = len(S)
p = ring_palette_bound(cls, n * k, S, k)
q = -(-m // p)
print(f"S: {m} points, palette {p} of {k}, forced class >= {q}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
dm, D2 = b.dim, b.D * b.D
bydist = defaultdict(list)
edges_in = []
for i in range(m):
    for j in range(i + 1, m):
        u, v = S[i], S[j]
        d = r[u] - r[v]
        sq = (b._field_square(d[None, :dm]) + b._field_square(d[None, dm:]))[0]
        if (min(u, v), max(u, v)) in Eset:
            edges_in.append((i, j))
        if not any(sq[1:]):
            bydist[Fr(int(sq[0]), D2)].append((i, j))
rat = sorted(bydist)
print(f"rational distances inside S: {[str(x) for x in rat]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def escapes(keys):
    pool = IDPool(start_from=m + 1)
    y = list(range(1, m + 1))
    f = list(CardEnc.atleast(lits=y, bound=q, vpool=pool,
                             encoding=EncType.seqcounter))
    for i, j in edges_in:
        f.append([-(i + 1), -(j + 1)])
    for key in keys:
        for i, j in bydist[key]:
            f.append([-(i + 1), -(j + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    out = s.solve()
    s.delete()
    return out


assert not escapes(rat), "the full list does not block -- nothing to minimise"
keep = list(rat)
for key in sorted(rat, key=lambda x: len(bydist[x])):
    cand = [x for x in keep if x != key]
    if cand and not escapes(cand):
        keep = cand
print(f"\nminimal blocking set: {len(keep)} distances "
      f"{[str(x) for x in keep]}  [{time.time()-t0:.0f}s]", flush=True)
for key in keep:
    print(f"   D={key} (d={float(key)**.5:.5f}), {len(bydist[key])} pairs, "
          f"closable: {bool(closable_distance(key))}", flush=True)
if len(keep) == 1:
    print("\n*** SINGLE DISTANCE: in every 5-colouring of G some pair at "
          f"distance^2 {keep[0]} inside S is monochromatic -- a FORCED PAIR",
          flush=True)
else:
    print(f"\na disjunction over {len(keep)} distances, not yet a forced pair",
          flush=True)
print("DONE", flush=True)
