"""Test the universe, not the seed: the cap is monotone, so one call decides
a whole space.

Every seed search here asked whether a particular subgraph has a capped ring,
and the space of seeds is 2^39 at one level and worse above.  That was the
wrong shape of question, because the property is MONOTONE.

  - "ring R cannot show all k colours" is preserved by ADDING vertices: more
    vertices means fewer colourings, so a ring that could not show them still
    cannot.
  - Conversely if the full point set's ring CAN show all k colours, then any
    subgraph containing that same ring can too -- restrict the colouring.

So for any ring of at least k points, testing the maximal point set answers the
question for EVERY subgraph that contains that ring.  A universe that fails
rules out its entire subgraph lattice at once, which is 2^n answers for the
price of one.

(Subgraphs whose ring is smaller are not covered, but a ring with fewer than k
points is capped trivially and carries nothing.)

Built as richly as the arithmetic allows: every point a short unit walk reaches
using every unit direction the field offers, inside a disc, then closed
dihedrally so the rings about the origin are as populated as they can be.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, deque
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
RADIUS = float(sys.argv[2]) if len(sys.argv) > 2 else 2.6
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 9000
NROT = int(sys.argv[4]) if len(sys.argv) > 4 else 5
DEPTH = int(sys.argv[5]) if len(sys.argv) > 5 else 3
t0 = time.time()
K = Field((2, 3, 5, 7, 11))
half = K.rational(Fr(1, 2))
rot60 = Rotation(half, K.sqrt(3) * half)
ZERO = Point(K.zero(), K.zero())

# every unit direction the field gives cheaply: the rotations joining a ring
# turn one unit vector into a family, and each is a step the walk can take
dirs, seen_d = [], set()
base = Point(K.rational(1), K.zero())
Ds = [Fr(x, y) for y in range(1, 7) for x in range(1, 7 * y + 1)]
rots = [Rotation(half, K.sqrt(3) * half)]
for D in sorted(set(Ds)):
    if D == 1 or not closable_distance(D):
        continue
    v = K.rational(4) * K.rational(D) - K.rational(1)
    c = v.c
    if any(x for x in c[1:]):
        continue
    q = c[0]
    s = None
    for r in K._prod:
        t = q / r
        n2, d2 = t.numerator, t.denominator
        rn, rd = int(round(n2 ** .5)), int(round(d2 ** .5))
        if rn * rn == n2 and rd * rd == d2:
            s = K.rational(Fr(rn, rd)) if r == 1 else K.sqrt(r) * K.rational(Fr(rn, rd))
            break
    if s is None:
        continue
    rots.append(Rotation(K.rational(1 - Fr(1, 2) / D),
                         s * K.rational(Fr(1, 2) / D)))
    if len(rots) >= NROT:
        break
u = base
for _ in range(6):
    for R in rots:
        w = R(u)
        if w not in seen_d:
            seen_d.add(w)
            dirs.append(w)
    u = rot60(u)
print(f"{len(dirs)} unit directions from {len(rots)} rotations"
      f"  [{time.time()-t0:.0f}s]", flush=True)

pool, seen, fr = [ZERO], {ZERO}, deque([(ZERO, 0)])
while fr and len(pool) < CAP:
    p, depth = fr.popleft()
    if depth >= DEPTH:
        continue
    for v in dirs:
        q = Point(p.x + v.x, p.y + v.y)
        if q in seen:
            continue
        if float(q.x) ** 2 + float(q.y) ** 2 > RADIUS * RADIUS:
            continue
        seen.add(q)
        pool.append(q)
        fr.append((q, depth + 1))
print(f"universe: {len(pool)} points within radius {RADIUS}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

b = IntBasis.covering(pool)
r = b.rows(pool)
head = b.overflow_headroom(r)
print(f"headroom {head:.2e}  [{time.time()-t0:.0f}s]", flush=True)
assert head < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(pool)
print(f"{n} points, {len(E)} edges ({len(E)/n:.2f}/v)"
      f"  [{time.time()-t0:.0f}s]", flush=True)
dm, D2 = b.dim, b.D * b.D
zi = pool.index(ZERO)
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
grp = defaultdict(list)
for off in np.nonzero(ok)[0]:
    v = Fr(int(sq[off, 0]), D2)
    if v and v != 1 and closable_distance(v):
        grp[v].append(int(off))
rings = [(D, m) for D, m in sorted(grp.items(), key=lambda t: -len(t[1]))
         if len(m) >= k]
print(f"{len(rings)} rings of >= {k} points, biggest "
      f"{[(str(D), len(m)) for D, m in rings[:6]]}  [{time.time()-t0:.0f}s]",
      flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
s0 = Solver(name="cd15", bootstrap_with=cls)
base_ok = s0.solve()
s0.delete()
print(f"{k}-colourable: {base_ok}  [{time.time()-t0:.0f}s]", flush=True)
if not base_ok:
    print(f"*** {n} points NOT {k}-COLOURABLE -- chi >= {k+1} ***", flush=True)
    sys.exit()
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
capped = []
for i, (D, mem) in enumerate(rings):
    if not sv.solve(assumptions=[sel0 + i]):
        capped.append((str(D), len(mem)))
        print(f"*** CAPPED: D={D}, {len(mem)} points ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if i % 20 == 19:
        print(f"   {i+1}/{len(rings)} rings, {len(capped)} capped"
              f"  [{time.time()-t0:.0f}s]", flush=True)
sv.delete()
print(f"\n{len(capped)} of {len(rings)} rings capped below {k} -- and by "
      f"monotonicity this answers for EVERY subgraph containing the ring"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
