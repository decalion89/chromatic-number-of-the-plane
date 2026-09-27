"""Apply the census to unions of G, where it has never been applied.

Every negative about G so far was read off a forced pair, a gateway or a
palette -- all of them projections.  The census is the quantity the bite
actually consumes, and it has a gradient the others lack: Sa scores 10 of 715
at four colours, the bite tightens it to SEVEN, and everything at five colours
so far scores 855 of 855, dead flat.

G is 5-chromatic, so every rotated copy of it must use all five colours in any
colouring of a union, and copies that share vertices constrain each other.
That is the only mechanism here that has ever produced rigidity -- it is what
the bite is -- and it has been measured on G by forced pairs, which found
nothing, but never by the census, which is finer.

Rotations are ranked by contact first, since the bite's six edges were a
structural choice and the sweep on Sa found rotations with eighteen times
more.  Then the union is built cumulatively and censused at every step, so a
trend shows even when no step closes.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
P0 = build_G(K, as_graph=False)
b0 = IntBasis.covering(P0)
r0 = b0.rows(P0)
E0 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b0, r0))
deg = defaultdict(int)
for a, c in E0:
    deg[a] += 1
    deg[c] += 1
hub = max(range(len(P0)), key=lambda i: deg[i])
HUB = P0[hub]
print(f"G: {len(P0)} pts, {len(E0)} edges, hub {hub} degree {deg[hub]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

Ds = []
for den in range(1, 19):
    for num in range(1, 25 * den + 1):
        D = Fr(num, den)
        if D != 1 and D not in Ds and closable_distance(D):
            Ds.append(D)
print(f"{len(Ds)} rings to try  [{time.time()-t0:.0f}s]", flush=True)
ranked = []
for D in Ds:
    rot = rotation_joining(D, K)
    for sgn in (+1, -1):
        rr = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(HUB)
        seen = set(P0)
        img = [rr(p) for p in P0]
        U = list(P0) + [q for q in img if q not in seen]
        b = IntBasis.covering(U)
        r = b.rows(U)
        if b.overflow_headroom(r) >= 1.0:
            continue
        E = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
        ranked.append((len(E) - 2 * len(E0), len(U), D, sgn))
ranked.sort(reverse=True)
print(f"top rotations about the hub (cross edges, points, ring, dir):",
      flush=True)
for row in ranked[:8]:
    print(f"   {row}", flush=True)


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


PARTS = [p for p in partitions(list(range(7))) if len(p) <= k]


def census(P, label):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        print(f"{label}: headroom exceeded", flush=True)
        return None
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    ci = P.index(HUB)
    key = {tuple(r[i]): i for i in range(n)}
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), d2)
        if v and v != 1 and closable_distance(v):
            grp[v].append(int(off))
    pick = None
    for D in sorted(grp, key=lambda D: (D != 4, -len(grp[D]))):
        mem, ms = grp[D], set(grp[D])
        pairs = [(j, key[tuple(2 * r[ci] - r[j])]) for j in mem
                 if tuple(2 * r[ci] - r[j]) in key
                 and key[tuple(2 * r[ci] - r[j])] in ms
                 and key[tuple(2 * r[ci] - r[j])] > j]
        if len(pairs) >= 3:
            six = [x for pr in pairs[:3] for x in pr]
            six.sort(key=lambda i: math.atan2(float(P[i].y - P[ci].y),
                                              float(P[i].x - P[ci].x)))
            pick = (D, len(mem), [ci] + six)
            break
    if pick is None:
        print(f"{label}: no usable ring", flush=True)
        return None
    D, size, W = pick
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"*** {label}: {n} pts NOT {k}-COLOURABLE ***", flush=True)
        sv.delete()
        return "WIN"
    surv = sum(1 for part in PARTS
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(part)
                                        for e in blk]))
    sv.delete()
    print(f"{label}: {n} pts, {len(E)} edges, ring D={D} ({size}), "
          f"census {surv} of {len(PARTS)}  [{time.time()-t0:.0f}s]", flush=True)
    return surv


census(list(P0), "G alone")
U, seen = list(P0), set(P0)
for cross, npts, D, sgn in ranked[:6]:
    rot = rotation_joining(D, K)
    rr = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(HUB)
    fresh = [q for q in (rr(p) for p in P0) if q not in seen]
    seen.update(fresh)
    U.extend(fresh)
    res = census(U, f"+ D={D}{'+' if sgn > 0 else '-'} ({cross} cross)")
    if res == "WIN" or res is None:
        break
print("DONE", flush=True)
