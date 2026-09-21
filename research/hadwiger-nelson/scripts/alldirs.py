"""How many unit directions does the group actually have?

The withdrawn proof failed because the group contains unit vectors I never put
in: it is dense in the plane, so norm-one elements turn up as longer words and
as sums that land on the circle by accident.  That was the error.  It is also
the opportunity -- every one of them is an edge direction I was not using, and
density is what the whole search wants.

So count them.  Take the point set, find every pair at distance exactly one
with the exact test, and read off the distinct difference vectors.  If that
set is much larger than the generators I supplied, then rebuilding the walk
with ALL of them gives a denser graph than anything built here, and it costs
nothing but the enumeration.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import deque, Counter
from hn.field import Field
from hn.geometry import Point
from hn.fast import IntBasis, fast_edges_complete
from quotient import units

t0 = time.time()
K = Field((3, 5, 7, 11))
ZERO = Point(K.zero(), K.zero())


def walk_rows(U, R, cap):
    basis = IntBasis.covering(U + [ZERO])
    urows = basis.rows(U)
    dim, D = basis.dim, basis.D
    steps = np.concatenate([urows, -urows], axis=0)
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
    return basis, np.array(out, dtype=np.int64)


for rots, emax, R, cap in (((3, 4), 1, 2.0, 3000), ((3, 4), 1, 2.0, 7000),
                           ((4,), 2, 2.0, 7000), ((3, 4, 7), 1, 2.0, 7000)):
    U = units(rots, emax)
    basis, rows = walk_rows(U, R, cap)
    E = fast_edges_complete(basis, rows)
    dirs = set()
    for a, b in E:
        d = rows[b] - rows[a]
        dirs.add(min(d.tobytes(), (-d).tobytes()))
    print(f"\nrotations {rots} e<={emax}, R={R}: {len(U)} generators supplied, "
          f"{len(rows)} points, {len(E)} exact unit edges, "
          f"{len(E)/len(rows):.2f} per vertex  [{time.time()-t0:.0f}s]",
          flush=True)
    print(f"   distinct unit DIRECTIONS in the graph: {len(dirs)} "
          f"(supplied {len(U)}) -> {len(dirs)/max(len(U),1):.1f} times as many"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("\nDONE", flush=True)
