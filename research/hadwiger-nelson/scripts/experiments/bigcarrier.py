"""A carrier at the scale five colours would need, and the cap question.

The construction needs a ring CAPPED below k, and nothing at five colours has
ever been capped here -- but nothing at five colours has been tried at the
scale the ratios predict either.  Sa is 397 points for four; the same ratio on
a 509-vertex gadget wants tens of thousands, which is out of reach, but a few
thousand is not and has not been looked at.

Built the way Sa is, one level up: take the graph the rhombus growth produces,
close it under the same twelve-element dihedral group about the origin, and ask
whether any ring about that origin is capped below five.  The closure is about
six times the input because the growth is greedy and not symmetric, so a
thousand-point input gives a six-thousand-point carrier.

The screen is one SAT call per ring -- can the ring show all five colours at
once -- and the satisfiable answer is the fast one, so an uncapped sweep is
cheap.  Unit rings are excluded: their cap is the colourability of their
centre and says nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = 5
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
t0 = time.time()
rng = random.Random(3)
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))
rot60 = Rotation(half, K.sqrt(3) * half)
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


def candidates(P, tries):
    out, have = [], set(P)
    n = len(P)
    for _ in range(tries):
        A, B = P[rng.randrange(n)], P[rng.randrange(n)]
        if A == B:
            continue
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
            if q not in have:
                have.add(q)
                out.append(q)
    return out


P = build_Sa(K)
for step in range(STEPS):
    cand = candidates(P, 30000)
    if not cand:
        break
    gb = IntBasis.covering(P + cand)
    dim, D2 = gb.dim, gb.D * gb.D
    Pr = gb.rows(P)
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(gb, Pr)))
    nb = [set() for _ in range(len(P))]
    for a, c in E:
        nb[a].add(c)
        nb[c].add(a)
    scored = []
    for q in cand:
        o = gb.rows([q])[0]
        dd = Pr - o
        s = gb._field_square(dd[:, :dim]) + gb._field_square(dd[:, dim:])
        hit = s[:, 0] == D2
        for m in range(1, dim):
            hit &= s[:, m] == 0
        nbq = set(int(x) for x in np.nonzero(hit)[0])
        if len(nbq) < 2:
            continue
        gain = 0
        nl = sorted(nbq)
        for ii in range(len(nl) - 1):
            di = Pr[nl[ii + 1:]] - Pr[nl[ii]]
            s3 = gb._field_square(di[:, :dim]) + gb._field_square(di[:, dim:])
            g3 = s3[:, 0] == 3 * D2
            for m in range(1, dim):
                g3 &= s3[:, m] == 0
            for off in np.nonzero(g3)[0]:
                gain += len(nb[nl[ii]] & nb[nl[ii + 1 + int(off)]])
        scored.append((gain, len(nbq), q))
    scored.sort(key=lambda t: (-t[0], -t[1]))
    P = list(dict.fromkeys(P + [q for g, dg, q in scored[:40] if g > 0]))
print(f"grown to {len(P)} points  [{time.time()-t0:.0f}s]", flush=True)

seen, U = set(), []
for p in P:
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                U.append(q)
            q = rot60(q)
print(f"dihedral closure: {len(U)} points  [{time.time()-t0:.0f}s]", flush=True)

b = IntBasis.covering(U)
r = b.rows(U)
head = b.overflow_headroom(r)
assert head < 1.0, f"headroom {head}"
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(U)
print(f"{n} points, {len(E)} edges ({len(E)/n:.2f}/v), headroom {head:.2e}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
dm, D2 = b.dim, b.D * b.D
zi = U.index(ZERO)
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
         if len(m) >= 4]
print(f"{len(rings)} rings about the origin, biggest "
      f"{[(str(D), len(m)) for D, m in rings[:6]]}  [{time.time()-t0:.0f}s]",
      flush=True)

cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
nv = n * k
s0 = Solver(name="cd15", bootstrap_with=cls)
base_ok = s0.solve()
s0.delete()
print(f"{k}-colourable: {base_ok}  [{time.time()-t0:.0f}s]", flush=True)
if not base_ok:
    print(f"*** {n} points NOT {k}-COLOURABLE -- chi >= {k+1} ***", flush=True)
    sys.exit()
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
        print(f"*** CAPPED: ring D={D}, {len(mem)} points ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
sv.delete()
print(f"\n{len(capped)} of {len(rings)} rings capped below {k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
