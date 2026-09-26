"""Turn a capped set into a forced pair at a named distance.

G carries a set S whose palette is 4 out of 5: in every proper 5-colouring, S
shows at most four colours.  Its 51 points are therefore split among at most
four colour classes, so one class meets S in at least ceil(|S|/4) points, and
a colour class is an independent set.

That is a statement about SOME subset of S, not a named one, and on its own it
names no distance.  But it can be made to name one.  If every independent
subset of S of that size contains a pair at distance exactly d, then whichever
class the pigeonhole produces, it contains such a pair -- and a monochromatic
pair at a known distance is exactly what a spindle consumes.

So the question is combinatorial and finite: is there an independent subset of
S, of the forced size, containing no pair at distance d?  Satisfiable means
the argument fails for that d.  Unsatisfiable means a forced pair at distance
d in every 5-colouring of G, which spindles to a sixth colour whenever the
field holds the joining rotation.

Asked for every distance that occurs inside S.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
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
S, _pts = pickle.load(open(SC + f"grown_{k}_{TARGET}.pkl", "rb"))
S = sorted(S)
m = len(S)
p = ring_palette_bound(cls, n * k, S, k)
q = -(-m // p)
print(f"S: {m} points, palette {p} of {k}, so some colour class meets S in "
      f"at least {q} points  [{time.time()-t0:.0f}s]", flush=True)
assert p <= TARGET

dm, D2 = b.dim, b.D * b.D
dist = {}
bydist = defaultdict(list)
for i in range(m):
    for j in range(i + 1, m):
        u, v = S[i], S[j]
        d = r[u] - r[v]
        sq = (b._field_square(d[None, :dm]) + b._field_square(d[None, dm:]))[0]
        key = (Fr(int(sq[0]), D2) if not any(sq[1:])
               else ("irr", tuple(int(x) for x in sq)))
        dist[(i, j)] = key
        bydist[key].append((i, j))
edges_in = [(i, j) for (i, j), key in dist.items()
            if (min(S[i], S[j]), max(S[i], S[j])) in Eset]
print(f"{len(bydist)} distinct distances inside S, {len(edges_in)} of the "
      f"{m*(m-1)//2} pairs are edges  [{time.time()-t0:.0f}s]", flush=True)
ranked = sorted(bydist.items(), key=lambda t: -len(t[1]))
print("the commonest distances and how many pairs carry them:", flush=True)
for key, prs in ranked[:10]:
    lab = str(key) if not isinstance(key, tuple) else "irrational"
    print(f"   {lab}: {len(prs)} pairs", flush=True)


def escapes(avoid_pairs):
    """Is there an independent q-subset of S avoiding all those pairs?"""
    pool = IDPool(start_from=m + 1)
    y = list(range(1, m + 1))
    f = list(CardEnc.atleast(lits=y, bound=q, vpool=pool,
                             encoding=EncType.seqcounter))
    for i, j in edges_in:                      # a colour class is independent
        f.append([-(i + 1), -(j + 1)])
    for i, j in avoid_pairs:
        f.append([-(i + 1), -(j + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    out = s.solve()
    s.delete()
    return out


print(f"\nasking, for each distance, whether an independent {q}-subset can "
      f"avoid it entirely:", flush=True)
wins = []
for key, prs in ranked:
    if isinstance(key, tuple):
        continue
    if not escapes(prs):
        wins.append((key, len(prs)))
        cl = closable_distance(key)
        print(f"*** NO ESCAPE at distance^2 {key} ({len(prs)} pairs): every "
              f"forced class contains such a pair, so in EVERY 5-colouring "
              f"of G some pair at that distance is monochromatic."
              f"  closable: {bool(cl)}  [{time.time()-t0:.0f}s]", flush=True)
# and the union of all rational distances, as a disjunction
allrat = [pr for key, prs in bydist.items() if not isinstance(key, tuple)
          for pr in prs]
print(f"\nwith every rational distance forbidden at once "
      f"({len(allrat)} pairs): escape = {escapes(allrat)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"distances with no escape: {len(wins)}", flush=True)
print("DONE", flush=True)
