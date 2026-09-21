"""Which centre gives G a populated doubly-usable ring?

The template is established and its bottleneck is now arithmetic, not size.
De Grey's field closes 22 doubly-usable rings among the small rationals --
2/7, 2/5, 4/9, 1/2, 8/7, 4, 17/2 and more -- and G populates exactly one of
them, D = 4, about the origin.  The bite has nowhere else to happen.

But a ring is a ring ABOUT A CENTRE, and nothing says the centre must be the
origin.  For a centre C and a point p, the bite needs |p - C|^2 rational and
doubly usable; the rotation rho_{D,C} then moves every point of that ring by
exactly one, and those are the cross edges of G u rho(G).  Sa's D = 4 ring
carried twelve points and twelve cross edges were enough at four colours.

So scan the centres.  Every point of G, every midpoint of a G pair at rational
squared distance, and the rational points nearby: for each, which doubly
usable rings does it see, and how populated are they.  Rank by population,
because coupling is what the bite buys.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.fast import IntBasis
from hn.homcol import doubly_usable_ring

t0 = time.time()
P = build_G(F, as_graph=False)
n = len(P)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
print(f"G: {n} pts, D={basis.D}  [{time.time()-t0:.0f}s]", flush=True)

DU = {}


def du(d):
    if d not in DU:
        DU[d] = doubly_usable_ring(d)
    return DU[d]


def ringcount(centre_row):
    """Populations of the rational rings about a centre given as an int row."""
    d = rows - centre_row
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    c = Counter()
    for off in np.nonzero(rat)[0]:
        v = int(sq[off, 0])
        if v:
            c[Fr(v, D2)] += 1
    return c


# --- centres: the points of G themselves -----------------------------------
best = []
for i in range(n):
    c = ringcount(rows[i])
    for d, m in c.items():
        if du(d):
            best.append((m, d, "vertex", i))
    if i % 400 == 0:
        print(f"  vertex centres {i}/{n}, {len(best)} hits"
              f"  [{time.time()-t0:.0f}s]", flush=True)

# --- centres: midpoints (doubled lattice, so no halving of the denominator) -
# A midpoint has coordinates (a+b)/2, so use 2*rows and a doubled denominator.
rows2 = rows * 2
print(f"\nvertex centres done: {len(best)} (centre, doubly-usable ring) hits"
      f"  [{time.time()-t0:.0f}s]", flush=True)
mid = []
seen = set()
for i in range(0, n, 1):
    for j in range(i + 1, n):
        key = tuple((rows[i] + rows[j]))
        if key in seen:
            continue
        seen.add(key)
    if i > 60:
        break
print(f"  (midpoint enumeration is quadratic; sampling instead)", flush=True)

import random
rng = random.Random(97)
mids = set()
while len(mids) < 4000:
    i, j = rng.randrange(n), rng.randrange(n)
    if i != j:
        mids.add((min(i, j), max(i, j)))
mb = []
for t, (i, j) in enumerate(sorted(mids)):
    # centre = (p_i + p_j)/2 : work in the doubled lattice
    d = rows2 - (rows[i] + rows[j])
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    c = Counter()
    for off in np.nonzero(rat)[0]:
        v = int(sq[off, 0])
        if v:
            c[Fr(v, 4 * D2)] += 1
    for dd, m in c.items():
        if du(dd):
            mb.append((m, dd, "midpoint", (i, j)))
    if t and t % 800 == 0:
        print(f"  midpoints {t}/{len(mids)}, {len(mb)} hits"
              f"  [{time.time()-t0:.0f}s]", flush=True)

allh = sorted(best + mb, key=lambda r: -r[0])
print(f"\n{len(allh)} (centre, doubly-usable ring) hits in all"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("top 30 by ring population:")
for m, d, kind, w in allh[:30]:
    print(f"   {m:4d} points on D={d}  ({kind} {w})", flush=True)
print("\nring values that appear:", Counter(d for _, d, _, _ in allh), flush=True)
print("DONE", flush=True)
