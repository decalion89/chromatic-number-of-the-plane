"""Feed the discovered directions back in, and repeat.

The group has more unit directions than were supplied -- 87 against 54, 235
against 162 -- and those extra ones are edges the walk never used, because it
only ever stepped along the generators it was given.  So feed them back: walk,
read off every unit direction the exact edge test finds, walk again with all
of them, repeat until the set stops growing.

Each round the step set grows, so each point's neighbourhood fills out further
and the density rises.  The limit is the set of ALL unit vectors of the group
inside the region, which is what the unit-distance graph on that region really
has -- so this converges to the right object rather than to an artefact of the
generators.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import deque
from hn.field import Field
from hn.geometry import Point
from hn.fast import IntBasis, fast_edges_complete
from quotient import units
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
ZERO = Point(K.zero(), K.zero())


def walk(basis, steprows, R, cap):
    dim, D = basis.dim, basis.D
    steps = np.concatenate([steprows, -steprows], axis=0)
    lim = R * R * D * D

    def norms(V):
        return (basis._field_square(V[:, :dim])
                + basis._field_square(V[:, dim:])) @ basis.sqrts

    z = np.zeros(2 * dim, dtype=np.int64)
    seen = {z.tobytes()}
    out = [z]
    q = deque([0])
    while q and len(out) < cap:
        i = q.popleft()
        cand = steps + out[i]
        for w in cand[norms(cand) <= lim]:
            b = w.tobytes()
            if b not in seen and len(out) < cap:
                seen.add(b)
                out.append(w)
                q.append(len(out) - 1)
    return np.array(out, dtype=np.int64)


def solve(n, E):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


for rots, emax, R, cap in (((3, 4), 1, 2.0, 7000), ((3, 4, 7), 1, 2.0, 7000)):
    U = units(rots, emax)
    basis = IntBasis.covering(U + [ZERO])
    steprows = basis.rows(U)
    print(f"\nrotations {rots}, R={R}, cap={cap}: starting from {len(U)} "
          f"directions  [{time.time()-t0:.0f}s]", flush=True)
    for rnd in range(1, 6):
        rows = walk(basis, steprows, R, cap)
        E = fast_edges_complete(basis, rows)
        dirs = {}
        for a, b in E:
            d = rows[b] - rows[a]
            key = min(d.tobytes(), (-d).tobytes())
            if key not in dirs:
                dirs[key] = d
        n = len(rows)
        ok = solve(n, sorted(set((min(a, b), max(a, b)) for a, b in E)))
        print(f"   round {rnd}: {len(steprows)} steps in, {n} pts, {len(E)} "
              f"edges, {len(E)/n:.2f} per vertex, {len(dirs)} directions out "
              f"-> {'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if not ok:
            print("*** a subgraph of the plane needs six colours ***",
                  flush=True)
            sys.exit(0)
        new = np.array(list(dirs.values()), dtype=np.int64)
        if len(new) <= len(steprows):
            print("   directions stopped growing", flush=True)
            break
        steprows = new
print("\nDONE", flush=True)
