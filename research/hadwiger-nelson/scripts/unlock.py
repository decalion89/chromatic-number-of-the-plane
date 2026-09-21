"""Which radical unlocks G's richest ring?

The template's bottleneck is arithmetic: the field closes twenty-two
doubly-usable rings and G populates one of them, D = 4, with twelve points --
Sa's ring, inherited and never improved on.  Thickening it by iterating the
bite adds twelve points and 2372 vertices per step, which is expensive.

Adjoining a radical is free.  The POINTS do not move: G is what it is, and
only the biting rotation rho_D needs sqrt(4D - 1) and the spindle
sqrt(16D - 1).  Adjoin the radical those want and a ring G already carries
becomes usable, at no cost in vertices at all.

So census every ring, not just the usable ones.  For each vertex of G as
centre, every rational ring and how many points sit on it; for each, the two
radicands the template demands; and then, over candidate adjunctions, which
one turns the biggest ring G HAS into a ring the field can bite.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.fast import IntBasis
from hn.homcol import closing_radicand

t0 = time.time()
P = build_G(F, as_graph=False)
n = len(P)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
print(f"G: {n} pts  [{time.time()-t0:.0f}s]", flush=True)

# best population of each rational ring value, over all vertex centres
best = defaultdict(int)
where = {}
for i in range(n):
    d = rows - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    c = Counter()
    for off in np.nonzero(rat)[0]:
        v = int(sq[off, 0])
        if v:
            c[Fr(v, D2)] += 1
    for dd, m in c.items():
        if m > best[dd]:
            best[dd] = m
            where[dd] = i
    if i % 400 == 0:
        print(f"  centre {i}/{n}, {len(best)} ring values"
              f"  [{time.time()-t0:.0f}s]", flush=True)

print(f"\n{len(best)} distinct rational ring values over all vertex centres"
      f"  [{time.time()-t0:.0f}s]", flush=True)
top = sorted(best.items(), key=lambda t: -t[1])[:28]
print("\nG's richest rings, and the two radicands the template demands:")
print(f"{'D':>10} {'pts':>5} {'centre':>7}  {'4D-1':>10} {'16D-1':>10}")
for d, m in top:
    a, b = closing_radicand(d), closing_radicand(4 * d)
    print(f"{str(d):>10} {m:5d} {where[d]:7d}  {'v'+str(a):>10} "
          f"{'v'+str(b):>10}", flush=True)


def reach(rads):
    out = {1}
    for p in rads:
        out |= {t * p for t in out}
    return out


BASE = (3, 5, 7, 11)
print("\nWhat one adjunction buys, ranked by the ring it unlocks:")
cands = [2, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
rows_out = []
for q in cands:
    R = reach(BASE + (q,))
    R0 = reach(BASE)
    gained = [(m, d) for d, m in best.items()
              if d != 1
              and closing_radicand(d) in R and closing_radicand(4 * d) in R
              and not (closing_radicand(d) in R0
                       and closing_radicand(4 * d) in R0)]
    gained.sort(reverse=True)
    rows_out.append((gained[0][0] if gained else 0, q, len(gained),
                     [(str(d), m) for m, d in gained[:5]]))
rows_out.sort(reverse=True)
for m, q, cnt, sample in rows_out:
    print(f"  adjoin sqrt({q:2d}): unlocks {cnt:3d} ring values, best {m:4d} "
          f"points  {sample}", flush=True)
print("\nDONE", flush=True)
