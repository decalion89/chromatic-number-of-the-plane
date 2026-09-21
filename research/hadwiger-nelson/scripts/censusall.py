"""The joint centre+ring census, as a comparable number across graphs.

De Grey's derivation runs on one number: how many patterns of a centre and six
ring points survive at k colours.  Sa scores 10 of 715 at four, which is
enough to eliminate everything, and 855 of 855 at five, which is nothing.  The
palette bound is a projection of this and throws away the shapes; the census
keeps them, and it is the quantity the bite actually consumes.

It is also comparable.  Six ring points and a centre is seven points whatever
the graph, so the count can be put side by side for every object built here
and the smallest one wins.  Six points are chosen antipodally closed -- three
pairs -- so the shapes mean the same thing they mean in Sa.

Cheap, too: the surviving calls are satisfiable and fast, and only the
eliminated ones cost, so a rigid graph is EXPENSIVE to census and a loose one
is instant.  The cost is itself a reading.
"""
import sys, time, math
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_S, build_Sa, build_Sb, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K, Point, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()
half = K.rational(Fr(1, 2))
rot30 = Rotation(K.sqrt(3) * half, half)
rot60 = Rotation(half, K.sqrt(3) * half)


def closure(seed, rot, order):
    seen, out = set(), []
    for p in seed:
        for base in (p, Point(p.x, -p.y)):
            q = base
            for _ in range(order):
                if q not in seen:
                    seen.add(q)
                    out.append(q)
                q = rot(q)
    return out


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


PARTS = [p for p in partitions(list(range(7)))]


def census(P, label, ks):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        print(f"{label}: headroom exceeded, skipped", flush=True)
        return
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    key = {tuple(r[i]): i for i in range(n)}
    best = None
    for ci in range(n):
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
        for D, mem in grp.items():
            ms = set(mem)
            pairs = []
            for j in mem:
                a = key.get(tuple(2 * r[ci] - r[j]))
                if a is not None and a in ms and a > j:
                    pairs.append((j, a))
            if len(pairs) >= 3:
                sixe = [x for pr in pairs[:3] for x in pr]
                sixe.sort(key=lambda i: math.atan2(float(P[i].y - P[ci].y),
                                                   float(P[i].x - P[ci].x)))
                cand = (ci, D, len(mem), [ci] + sixe)
                if best is None or D == 4:
                    best = cand
                if D == 4:
                    break
        if best and best[1] == 4:
            break
    if best is None:
        print(f"{label}: {n} pts, no ring with three antipodal pairs",
              flush=True)
        return
    ci, D, size, W = best
    for k in ks:
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, c in E:
            for col in range(k):
                cls.append([-(1 + a * k + col), -(1 + c * k + col)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        if not sv.solve():
            print(f"{label}: {n} pts *** NOT {k}-COLOURABLE ***", flush=True)
            sv.delete()
            continue
        tot = surv = alone = 0
        for part in PARTS:
            if len(part) > k:
                continue
            tot += 1
            lits = [1 + W[e] * k + bi
                    for bi, blk in enumerate(part) for e in blk]
            if sv.solve(assumptions=lits):
                surv += 1
                if any(blk == [0] for blk in part):
                    alone += 1
        sv.delete()
        print(f"{label}: {n} pts, centre {ci}, ring D={D} ({size} pts), "
              f"k={k}: {surv} of {tot} patterns, centre alone in {alone}"
              f"  [{time.time()-t0:.0f}s]", flush=True)


S = build_S(K)
census(build_Sa(K), "Sa (D6 closure)", (4, 5))
census(closure(S, rot30, 12), "D12 closure of S", (4, 5))
Y = build_Y(K)
census(Y, "Y (the bite)", (4, 5))
census(closure(Y, rot30, 12), "D12 closure of Y", (5,))
census(build_G(K, as_graph=False), "G", (5,))
print("DONE", flush=True)
