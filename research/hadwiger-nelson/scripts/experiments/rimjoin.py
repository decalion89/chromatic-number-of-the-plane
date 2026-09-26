"""Join two copies of G at the rim, where its fifth colour actually lives.

Two measurements point the same way.  Adding 7206 points inside a radius-1.0
disc left the chromatic number at 4, so filling in buys nothing.  And the disc
of radius 2.0 holds 1569 of G's 1581 points while still colouring with four,
so the fifth colour rests entirely on the twelve points outside it -- six
unit-distance pairs on the rim, degrees 6 and 7 against a mean of 9.96.

What a higher chromatic number needs here is reach, and the reach ends in tip
pairs.  So the move is to put two copies of G rim to rim and stitch them
there, rather than to thicken either one.

Candidate placements come from the graph itself, which keeps everything
exact: translate by the difference of two points, then by a unit vector, so
that a chosen rim point of one copy lands a unit from a chosen rim point of
the other.  Ranked by how many cross edges the placement creates, and by how
many of those land on rim points, since those are the constrained ones.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
SC_PATH = ("/tmp/hn/")
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver
sys.path.insert(0, SC_PATH)
from pairsat import pairs_at

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
TRIES = int(sys.argv[2]) if len(sys.argv) > 2 else 12
t0 = time.time()
SC = SC_PATH
G = build_G(K, as_graph=False)
n0 = len(G)
cx = sum(float(p.x) for p in G) / n0
cy = sum(float(p.y) for p in G) / n0
rim = [i for i, p in enumerate(G)
       if (float(p.x) - cx) ** 2 + (float(p.y) - cy) ** 2 > 4.0]
print(f"G: {n0} points, rim of {len(rim)}  [{time.time()-t0:.0f}s]",
      flush=True)
b0 = IntBasis.covering(G)
r0 = b0.rows(G)
assert b0.overflow_headroom(r0) < 1.0
E0 = [(a, c) for a, c in fast_edges_complete(b0, r0) if a < c]
# unit vectors available in the graph, as integer rows
units = []
seen = set()
for a, c in E0:
    for d in (r0[c] - r0[a], r0[a] - r0[c]):
        t = tuple(int(x) for x in d)
        if t not in seen:
            seen.add(t)
            units.append(np.array(d))
    if len(units) > 600:
        break
print(f"{len(units)} unit vectors  [{time.time()-t0:.0f}s]", flush=True)

rimrows = r0[rim]
cand = {}
for i in rim:
    for j in rim:
        base = r0[i] - r0[j]
        for e in units[::3]:
            t = tuple(int(x) for x in (base + e))
            cand[t] = cand.get(t, 0) + 1
print(f"{len(cand)} candidate translations  [{time.time()-t0:.0f}s]",
      flush=True)
best = sorted(cand.items(), key=lambda kv: -kv[1])[:TRIES]
rimset = set(rim)
for trow, score in best:
    tv = np.array(trow, dtype=np.int64)
    rows = np.vstack([r0, r0 + tv])
    uniq = {}
    for i, rw in enumerate(rows):
        uniq.setdefault(tuple(int(x) for x in rw), i)
    keep = sorted(uniq.values())
    A = rows[keep]
    hr = b0.overflow_headroom(A)
    if hr >= 1.0:
        continue
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b0, A)))
    n = len(A)
    orig = {i for i, v in enumerate(keep) if v < n0}
    cross = sum(1 for a, c in E if (a in orig) != (c in orig))
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    for col in range(1, k):
        cls.append([-(1 + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    # The union doubles the points and its interface is rich, so its
    # two-distance graph is far denser than either copy's -- and there the
    # target is six, not five.
    E2 = set(pairs_at(b0, A, int(Fr(3) * b0.D * b0.D)))
    both = sorted(set(E) | E2)
    ch2 = None
    for kk in range(5, 9):
        c2 = [[1 + v * kk + c for c in range(kk)] for v in range(n)]
        for a, c in both:
            for col in range(kk):
                c2.append([-(1 + a * kk + col), -(1 + c * kk + col)])
        for col in range(1, kk):
            c2.append([-(1 + col)])
        sv = Solver(name="cd15", bootstrap_with=c2)
        good = sv.solve()
        sv.delete()
        if good:
            ch2 = kk
            break
    mark = "   <<< SIX: forced pair at sqrt3" if (ch2 or 0) >= 6 else ""
    print(f"placement (rim score {score}): {n} points, {len(E)} edges, "
          f"{cross} cross -> "
          f"{'%d-colourable' % k if ok else '*** NOT %d-COLOURABLE ***' % k}"
          f"; two-distance chi {ch2}{mark}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        pickle.dump([Point(*b0.to_field(A[i])) if hasattr(b0, "to_field")
                     else None for i in range(n)],
                    open(SC + "WITNESS_rim.pkl", "wb"))
        print("   witness rows saved", flush=True)
        np.save(SC + "WITNESS_rim_rows.npy", A)
        break
print("DONE", flush=True)
