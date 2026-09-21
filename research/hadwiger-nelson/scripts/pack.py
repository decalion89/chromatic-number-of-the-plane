"""Pack five-colour gadgets by overlapping G with its own translates.

The account of why five colours stalls is gadget density.  Sa is saturated
with Moser spindles -- 576 of them on 397 points, 1.45 per point -- so its
4-colourings are pinned; G contains essentially ONE 5-critical subgraph, about
0.001 per point, and its 5-colourings are free in every direction measured.

Density can be bought directly.  A translate G + t contains a copy of G, and
if t is chosen so the two overlap heavily, the union is barely larger than G
while containing two 5-chromatic subgraphs that share most of their vertices.
Each must use all five colours, and they constrain each other through the
shared part.  Iterating over several such translations packs many copies into
little room, which is exactly what Sa does with spindles.

The translations worth using are the frequent differences of G's own points:
one pass over all 1581^2 differences, counted by hashing, gives the overlap of
every candidate at once.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
P = build_G(K, as_graph=False)
n = len(P)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E0 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
print(f"G: {n} pts, {len(E0)} edges  [{time.time()-t0:.0f}s]", flush=True)

# every difference, counted -- the count IS the overlap of that translate
cnt = Counter()
for i in range(n):
    d = r - r[i]
    for row in d:
        cnt[row.tobytes()] += 1
zero = np.zeros(2 * b.dim, dtype=np.int64).tobytes()
del cnt[zero]
top = cnt.most_common(40)
print(f"{len(cnt)} distinct translations; best overlaps "
      f"{[c for _, c in top[:12]]}  [{time.time()-t0:.0f}s]", flush=True)


def as_point(buf):
    v = np.frombuffer(buf, dtype=np.int64)
    inv = Fr(1, b.D)
    x = K.element([Fr(int(v[m])) * inv for m in range(b.dim)])
    y = K.element([Fr(int(v[b.dim + m])) * inv for m in range(b.dim)])
    return Point(x, y)


def study(U, label):
    bb = IntBasis.covering(U)
    rr = bb.rows(U)
    if bb.overflow_headroom(rr) >= 1.0:
        print(f"  {label}: headroom exceeded", flush=True)
        return None
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(bb, rr)))
    m = len(U)
    cls = [[1 + v * k + c for c in range(k)] for v in range(m)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    print(f"  {label}: {m} pts, {len(E)} edges ({len(E)/m:.2f}/v), "
          f"{k}-colourable={ok}  [{time.time()-t0:.0f}s]", flush=True)
    return ok


U, seen = list(P), set(P)
used = 0
for buf, overlap in top:
    t = as_point(buf)
    fresh = [q for q in (Point(p.x + t.x, p.y + t.y) for p in P)
             if q not in seen]
    seen.update(fresh)
    U.extend(fresh)
    used += 1
    ok = study(U, f"+ translate {used} (overlap {overlap}, {len(fresh)} new)")
    if ok is None:
        break
    if not ok:
        print(f"*** {len(U)} points NOT {k}-COLOURABLE -- chi >= {k+1} ***",
              flush=True)
        break
    if len(U) > 14000:
        print("  size cap reached", flush=True)
        break
print("DONE", flush=True)
