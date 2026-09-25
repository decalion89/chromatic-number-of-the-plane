"""Edges straight from the generators -- a subgraph is enough.

Density is collected in a second pass over every point, not
only the ones the queue reached; see the note in walk().

Confining twenty thousand points to a disc of radius 1.5 defeats the cell
bucketing that fast_edges_complete relies on: cells have side one, so nearly
every point lands in the same few cells and the search degenerates to all
pairs.  It stalled there.

It is also unnecessary.  By construction p and p + u are at distance exactly
one for every generator u, so those edges come free from the walk itself, in
O(n |U|) rather than O(n^2).  What that builds is a SUBGRAPH of the true
unit-distance graph -- it can miss coincidental unit pairs -- and a subgraph is
all the argument needs: if the subgraph is not 5-colourable then neither is the
graph, and neither is the plane.  A subgraph that colours proves nothing about
the full one, which is stated rather than hidden.

Density is now a dial: every interior point has |±U| neighbours, so the graph
approaches |U| edges per vertex, against 4.98 for the whole de Grey family and
6.59 for the best sumset.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import deque
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from hn.fast import IntBasis
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())
W = _rot60(K)
ROT = {D: rotation_joining(Fr(D), K) for D in (2, 3, 4, 7)}


def units(rots, emax):
    out, seen = [], set()
    for a in range(6):
        for es in itertools.product(*[range(-emax, emax + 1) for _ in rots]):
            p = ONE
            for _ in range(a):
                p = W(p)
            for D, e in zip(rots, es):
                r = ROT[D]
                use = r if e >= 0 else Rotation(r.cos, -r.sin)
                for _ in range(abs(e)):
                    p = use(p)
            if p not in seen:
                seen.add(p)
                out.append(p)
    return out


def walk(U, R, cap):
    """BFS inside radius R, keeping the generator edges as we go."""
    basis = IntBasis.covering(U + [Point(K.zero(), K.zero())])
    rows = basis.rows(U)
    dim, D = basis.dim, basis.D
    steps = np.concatenate([rows, -rows], axis=0)
    lim = R * R * D * D

    def norms(V):
        return (basis._field_square(V[:, :dim])
                + basis._field_square(V[:, dim:])) @ basis.sqrts

    z = np.zeros(2 * dim, dtype=np.int64)
    idx = {z.tobytes(): 0}
    out = [z]
    q = deque([0])
    E = set()
    while q and len(out) < cap:
        i = q.popleft()
        cand = steps + out[i]
        keep = cand[norms(cand) <= lim]
        for w in keep:
            b = w.tobytes()
            j = idx.get(b)
            if j is None:
                if len(out) >= cap:
                    continue
                j = len(out)
                idx[b] = j
                out.append(w)
                q.append(j)
            if i != j:
                E.add((min(i, j), max(i, j)))
    # A second pass over EVERY point, not just the ones the queue reached.
    # Collecting edges only while processing meant the points added last --
    # the whole frontier, and with many generators that is most of them --
    # kept none of their own.  That is why density FELL as generators were
    # added: 18 generators gave 4.66 per vertex and 162 gave 1.75, which is
    # backwards.  This pass is O(n |U|) and fixes it.
    arr = np.array(out, dtype=np.int64)
    for u in steps:
        for i, w in enumerate(arr + u):
            j = idx.get(w.tobytes())
            if j is not None and i != j:
                E.add((min(i, j), max(i, j)))
    return len(out), sorted(E)


def solve(n, E):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


# Density rises with the CAP, not with the generator count: at a fixed radius
# more points means more of them have their whole neighbourhood present, and
# the limit is |U| edges per vertex.  With thirty generators that limit is
# thirty, six times the de Grey family's 4.98.  So push the cap, on the two
# configurations that did best.
for rots, emax in (((4,), 2), ((4,), 1)):
    U = units(rots, emax)
    print(f"\nrotations {rots}, exponents <= {emax}: {len(U)} unit vectors, "
          f"degree up to {2*len(U)}, density limit {len(U)}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for R, cap in ((2.2, 45000), (2.2, 90000), (3.0, 90000)):
        n, E = walk(U, R, cap)
        if n < 60:
            continue
        print(f"   R={R} cap={cap}: {n} pts, {len(E)} edges, "
              f"{len(E)/n:.2f} per vertex; solving"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        ok = solve(n, E)
        print(f"      -> {'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if not ok:
            print("*** a subgraph of the plane needs six colours ***",
                  flush=True)
            sys.exit(0)
print("\nDONE", flush=True)
