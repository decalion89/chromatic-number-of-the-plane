"""The universe built the right way: unit-circle intersections, not walks.

The walk universe failed its own calibration -- no capped ring even at four
colours, where one is known -- because a sum of unit vectors never lands where
de Grey's points are.  His sit at radii 0.168, 0.264, 0.292, 0.333, and they
arise as INTERSECTIONS of unit circles about points already present, which is
the same generator the growth here uses.

So build the universe that way and check the calibration before believing
anything: the closure of a small start under "add both points at distance one
from these two", iterated inside a disc.  Then ask at FOUR colours whether any
ring is capped.  If the answer is still zero the universe is still wrong and
the five-colour question has not been reached; if it is not zero, the
instrument has found the known case and its five-colour reading means
something.

Seeded with Sa itself, so the calibration cannot fail for lack of the answer
being present -- and the interesting question is then how far the closure goes
beyond it.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 12000
RADIUS = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
t0 = time.time()
rng = random.Random(5)
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
print(f"start: Sa, {len(P)} points  [{time.time()-t0:.0f}s]", flush=True)
for rnd in range(ROUNDS):
    fresh = []
    n = len(P)
    order = list(range(n))
    rng.shuffle(order)
    for ii in order:
        A = P[ii]
        for jj in order:
            if jj <= ii:
                continue
            B = P[jj]
            D = A.dist2(B)
            f = float(D)
            if not .05 < f < 3.99:
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
            if len(P) + len(fresh) > CAP:
                break
        if len(P) + len(fresh) > CAP:
            break
    P = P + fresh
    print(f"round {rnd+1}: +{len(fresh)} -> {len(P)} points"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if len(P) >= CAP:
        break

b = IntBasis.covering(P)
r = b.rows(P)
head = b.overflow_headroom(r)
print(f"headroom {head:.2e}  [{time.time()-t0:.0f}s]", flush=True)
assert head < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{n} points, {len(E)} edges ({len(E)/n:.2f}/v)"
      f"  [{time.time()-t0:.0f}s]", flush=True)
dm, D2 = b.dim, b.D * b.D
zi = P.index(ZERO)
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
sv.delete()
print(f"\n{len(capped)} of {len(rings)} rings capped below {k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
