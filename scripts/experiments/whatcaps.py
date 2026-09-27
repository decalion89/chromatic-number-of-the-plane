"""What distinguishes a capped ring from an uncapped one?

Nobody has a criterion for which rings cap; de Grey found one ring and used it.
The calibrated universe changes that, because it produces POSITIVE AND NEGATIVE
examples inside a single graph: at four colours, four of its eleven rings are
capped and seven are not, all of them centred on the same vertex of the same
point set.  That is a classification problem with clean labels, and it has
never been posable before.

Size is already ruled out by inspection -- the largest ring, D = 3 on 42
points, is NOT capped, while D = 16/9 on six IS.  So the feature is structural.

Measured per ring: its population; the unit edges among its own points, which
make it a circulant whose step is the angle whose chord is one; how those
components look; how many antipodal pairs it carries; how many of its points
are adjacent to the centre; and the palette bound itself.  Then the two groups
are printed side by side and left to speak.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, doubly_usable_ring
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
RADIUS = 4.5
CAP = 12000
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))
ZERO = Point(K.zero(), K.zero())


def rsqrt(e):
    if any(x for x in e.c[1:]):
        return None
    v = e.c[0]
    if v <= 0:
        return None
    for rr in CLASSES:
        q = v / rr
        nn, dd = q.numerator, q.denominator
        rn, rd = int(round(nn ** .5)), int(round(dd ** .5))
        if rn * rn == nn and rd * rd == dd:
            s = K.rational(Fr(rn, rd))
            return s if rr == 1 else K.sqrt(rr) * s
    return None


P = build_Sa(K)
have = set(P)
fresh = []
n0 = len(P)
for i in range(n0):
    A = P[i]
    for j in range(i + 1, n0):
        B = P[j]
        D = A.dist2(B)
        if not .05 < float(D) < 3.99:
            continue
        sD, s4 = rsqrt(D), rsqrt(K.rational(4) - D)
        if sD is None or s4 is None:
            continue
        inv = K.rational(1) / sD
        mx, my = (A.x + B.x) * half, (A.y + B.y) * half
        nx = -(B.y - A.y) * inv * s4 * half
        ny = (B.x - A.x) * inv * s4 * half
        for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if q in have:
                continue
            if float(q.x) ** 2 + float(q.y) ** 2 > RADIUS * RADIUS:
                continue
            have.add(q)
            fresh.append(q)
P = P + fresh
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
Eset = set(E)
n = len(P)
dm, D2 = b.dim, b.D * b.D
zi = P.index(ZERO)
print(f"{n} points, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]", flush=True)
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
grp = defaultdict(list)
for off in np.nonzero(ok)[0]:
    v = Fr(int(sq[off, 0]), D2)
    if v and closable_distance(v):
        grp[v].append(int(off))
rings = [(D, m) for D, m in sorted(grp.items(), key=lambda t: -len(t[1]))
         if len(m) >= k]
key = {tuple(r[i]): i for i in range(n)}

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
nv = n * k
ind = [nv + 1 + c for c in range(k)]
enc = CardEnc.atleast(lits=list(ind), bound=k, top_id=nv + k,
                      encoding=EncType.seqcounter)
sel0 = max([nv + k] + [abs(x) for cl in enc.clauses for x in cl]) + 1
extra = list(enc.clauses)
for i, (D, mem) in enumerate(rings):
    for c in range(k):
        extra.append([-(sel0 + i), -ind[c]] + [1 + v * k + c for v in mem])
sv = Solver(name="cd15", bootstrap_with=cls + extra)
assert sv.solve(), "the universe must colour"
print(f"{'ring':>9} {'pts':>4} {'inner':>6} {'comps':>6} {'maxdeg':>7} "
      f"{'anti':>5} {'atcentre':>9} {'2xusable':>9}  capped", flush=True)
rows = []
for i, (D, mem) in enumerate(rings):
    ms = set(mem)
    inner = [(a, c) for a in mem for c in mem
             if a < c and (a, c) in Eset]
    nb = defaultdict(set)
    for a, c in inner:
        nb[a].add(c)
        nb[c].add(a)
    seen, comps = set(), 0
    for v in mem:
        if v in seen:
            continue
        comps += 1
        stack = [v]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(nb[x] - seen)
    anti = sum(1 for j in mem
               if key.get(tuple(2 * r[zi] - r[j])) in ms
               and key[tuple(2 * r[zi] - r[j])] > j)
    atc = sum(1 for j in mem if (min(zi, j), max(zi, j)) in Eset)
    capped = not sv.solve(assumptions=[sel0 + i])
    rows.append((D, len(mem), len(inner), comps,
                 max((len(nb[v]) for v in mem), default=0), anti, atc,
                 doubly_usable_ring(D), capped))
    print(f"{str(D):>9} {len(mem):>4} {len(inner):>6} {comps:>6} "
          f"{max((len(nb[v]) for v in mem), default=0):>7} {anti:>5} "
          f"{atc:>9} {str(doubly_usable_ring(D)):>9}  "
          f"{'YES' if capped else '.'}  [{time.time()-t0:.0f}s]", flush=True)
sv.delete()
cap = [x for x in rows if x[-1]]
unc = [x for x in rows if not x[-1]]
print(f"\ncapped {len(cap)}, uncapped {len(unc)}", flush=True)
for name, idx in (("points", 1), ("inner edges", 2), ("components", 3),
                  ("max inner degree", 4), ("antipodal pairs", 5),
                  ("adjacent to centre", 6)):
    if cap and unc:
        print(f"   {name:>20}: capped {sorted(x[idx] for x in cap)}  "
              f"uncapped {sorted(x[idx] for x in unc)}", flush=True)
print("DONE", flush=True)
