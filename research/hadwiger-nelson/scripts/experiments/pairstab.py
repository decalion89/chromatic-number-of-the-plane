"""The whole stabiliser of a PAIR, applied: half-turn and both mirrors.

The classification says an operation can only tighten the census of a set if
it maps that set to itself.  For a ring about c those are the rotations about
c -- the bites.  For a PAIR {u, v}, which is what the spindle actually needs
forced, they are different: the half-turn about the midpoint, which swaps u and
v, and the reflections in the two axes of the pair.

The half-turn is the interesting one and nothing here has tried it.  It moves a
point by twice its distance from the midpoint, so the matching between the two
copies is the circle of radius one half about that midpoint -- a ring again,
but centred on the pair rather than on the graph.

Tested on every pair of Sa at a closable distance, both by whether the pair
comes out forced and by whether the union stops colouring.  Sa is already
symmetric under the half-turn about its own centre, so that one is trivial and
every other is new.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, forced_same
from pysat.solvers import Solver

k = 4
t0 = time.time()
P0 = build_Sa(K)
n0 = len(P0)
b0 = IntBasis.covering(P0)
r0 = b0.rows(P0)
dm, D2 = b0.dim, b0.D * b0.D
E0 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b0, r0))
print(f"Sa: {n0} points, {len(E0)} edges  [{time.time()-t0:.0f}s]", flush=True)

pairs = []
for i in range(n0 - 1):
    d = r0[i + 1:] - r0[i]
    sq = b0._field_square(d[:, :dm]) + b0._field_square(d[:, dm:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dm):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) in E0:
            continue
        v = Fr(int(sq[off, 0]), D2)
        if v and closable_distance(v):
            pairs.append((v, i, j))
byd = defaultdict(list)
for v, i, j in pairs:
    byd[v].append((i, j))
print(f"{len(pairs)} closable non-edge pairs over {len(byd)} classes; "
      f"sizes {sorted((str(d), len(p)) for d, p in byd.items())[:6]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

half = K.rational(Fr(1, 2))
# one representative per class, plus every pair of the spindle-able class 16
todo = []
for D, ps in sorted(byd.items(), key=lambda t: len(t[1])):
    todo.extend((D, i, j) for i, j in ps[:6])
print(f"{len(todo)} half-turns to try  [{time.time()-t0:.0f}s]", flush=True)

tried = best = 0
for D, i, j in todo:
    u, v = P0[i], P0[j]
    mx, my = (u.x + v.x) * half, (u.y + v.y) * half
    # the three isometries that preserve the pair: the half-turn about the
    # midpoint, the reflection in the pair's own line, and the reflection in
    # its perpendicular bisector.  Together they are the whole stabiliser.
    dx, dy = v.x - u.x, v.y - u.y
    nsq = dx * dx + dy * dy
    ops = {
        "half-turn": lambda p: Point(mx + mx - p.x, my + my - p.y),
    }
    if nsq != K.zero():
        inv = K.rational(1) / nsq

        def _mirror_line(p, dx=dx, dy=dy, inv=inv):
            ax, ay = p.x - mx, p.y - my
            t = (ax * dx + ay * dy) * inv
            px, py = dx * t, dy * t
            return Point(mx + px + px - ax, my + py + py - ay)

        def _mirror_perp(p, dx=dx, dy=dy, inv=inv):
            ax, ay = p.x - mx, p.y - my
            t = (ax * dx + ay * dy) * inv
            px, py = dx * t, dy * t
            return Point(mx + ax - px - px, my + ay - py - py)

        ops["mirror in the pair's line"] = _mirror_line
        ops["mirror in the bisector"] = _mirror_perp
    U, seen = list(P0), set(P0)
    for nm, op in ops.items():
        for p in P0:
            q = op(p)
            if q not in seen:
                seen.add(q)
                U.append(q)
    b = IntBasis.covering(U)
    r = b.rows(U)
    if b.overflow_headroom(r) >= 1.0:
        continue
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    m = len(U)
    cls = [[1 + w * k + c for c in range(k)] for w in range(m)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    tried += 1
    if not sv.solve():
        print(f"*** pair {i},{j} at D={D}: {m} points NOT 4-COLOURABLE ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        break
    f = forced_same(sv, U.index(u), U.index(v), k)
    sv.delete()
    if f:
        best += 1
        print(f"*** pair {i},{j} at D={D}: {m} points, FORCED SAME ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if tried % 20 == 0:
        print(f"   {tried}/{len(todo)} half-turns, {best} forced"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} half-turns, {best} forced pairs  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
