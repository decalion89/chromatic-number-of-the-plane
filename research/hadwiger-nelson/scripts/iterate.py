"""Is the bite a trick or a process?

The census says the bite tightens: Sa scores 10 of 715 at four colours and Y,
which is Sa with one rotated copy, scores 7.  If biting AGAIN tightens again,
the bite is a repeatable operation and the chain has more levels in it than
de Grey used.  If it does not, the tightening is a one-off and the three
antipodal pairs are the end of that line.

Answering costs seven SAT calls, not seven hundred.  Adding vertices can only
REMOVE surviving patterns, never add one, so the survivors of a bigger graph
are a subset of the survivors of a smaller one.  Testing only Y's seven
survivors on the twice-bitten graph therefore settles it exactly.

Run down the chain: Y, then Y with one more bite, then two, as far as the
integer certificate allows.
"""
import sys, time, math
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Y
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
t0 = time.time()
rho = rotation_joining(4, K)
ZERO = Point(K.zero(), K.zero())


def partitions(seq):
    if not seq:
        yield []
        return
    first, rest = seq[0], seq[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p


ALL = [p for p in partitions(list(range(7))) if len(p) <= k]


def setup(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        return None
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    ci = P.index(ZERO)
    key = {tuple(r[i]): i for i in range(n)}
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
    ms = set(ring)
    pairs = [(j, key[tuple(2 * r[ci] - r[j])]) for j in ring
             if tuple(2 * r[ci] - r[j]) in key
             and key[tuple(2 * r[ci] - r[j])] in ms
             and key[tuple(2 * r[ci] - r[j])] > j]
    six = [x for pr in pairs[:3] for x in pr]
    six.sort(key=lambda i: math.atan2(float(P[i].y), float(P[i].x)))
    W = [ci] + six
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    return n, len(E), len(ring), W, cls


def survivors(state, cands):
    n, ne, rs, W, cls = state
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return "UNCOLOURABLE"
    out = []
    for part in cands:
        lits = [1 + W[e] * k + bi
                for bi, blk in enumerate(part) for e in blk]
        if sv.solve(assumptions=lits):
            out.append(part)
    sv.delete()
    return out


P = build_Y(K)
st = setup(P)
surv = survivors(st, ALL)
print(f"Y: {st[0]} pts, {st[1]} edges, ring {st[2]}: {len(surv)} of "
      f"{len(ALL)} patterns  [{time.time()-t0:.0f}s]", flush=True)
for p in surv:
    print(f"    {[sorted(b) for b in p]}", flush=True)

cur = list(P)
for level in range(1, 5):
    img = [rho(p) for p in cur]
    seen = set(cur)
    cur = cur + [q for q in img if q not in seen]
    st = setup(cur)
    if st is None:
        print(f"bite {level}: headroom exceeded, stop", flush=True)
        break
    nxt = survivors(st, surv)
    if nxt == "UNCOLOURABLE":
        print(f"*** bite {level}: {st[0]} pts NOT {k}-COLOURABLE ***",
              flush=True)
        break
    print(f"bite {level}: {st[0]} pts, {st[1]} edges, ring {st[2]}: "
          f"{len(nxt)} of {len(surv)} survive  [{time.time()-t0:.0f}s]",
          flush=True)
    for p in nxt:
        print(f"    {[sorted(b) for b in p]}", flush=True)
    if len(nxt) == len(surv):
        print("    no tightening -- the bite is a one-off", flush=True)
        break
    surv = nxt
print("DONE", flush=True)
