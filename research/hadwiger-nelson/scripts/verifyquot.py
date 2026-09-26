"""Check the periodic colouring on real points, not just on the quotient.

The quotient argument says: pick phi : Z^r -> Z_5 x Z_5 with no unit vector in
its kernel, 5-colour Cay(Z_5 x Z_5, phi(U)), and give every point of the group
the colour of its image.  Adjacent points differ by a unit vector, whose image
is non-zero and whose endpoints therefore differ in the quotient colouring.

That is a proof, and it is short enough to be wrong in a way that a SAT call
on the quotient would not catch -- a mistake in the coordinates, or in what
counts as a unit vector.  So take actual points of the group, build the actual
unit-distance graph on them with the exact edge test, colour them through phi,
and count monochromatic edges.  It must be zero.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, time, itertools, random
sys.path.insert(0, HN_DIR)
import numpy as np
from fractions import Fraction as Fr
from collections import deque
from hn.geometry import Point
from hn.field import Field
from hn.fast import IntBasis, fast_edges_complete
from quotient import units, coords
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
rng = random.Random(20240921)
A = (5, 5)


def find_phi(ints, r):
    for _ in range(200000):
        phi = [[rng.randrange(a) for a in A] for _ in range(r)]
        imgs = []
        ok = True
        for co in ints:
            g = tuple(sum(c * phi[t][s] for t, c in enumerate(co)) % A[s]
                      for s in range(len(A)))
            if not any(g):
                ok = False
                break
            imgs.append(g)
        if not ok:
            continue
        # 5-colour the quotient
        E = set()
        code = lambda v: v[0] * A[1] + v[1]
        for v in itertools.product(range(A[0]), range(A[1])):
            i = code(v)
            for g in set(imgs):
                j = code(((v[0] + g[0]) % A[0], (v[1] + g[1]) % A[1]))
                if i != j:
                    E.add((min(i, j), max(i, j)))
        size = A[0] * A[1]
        cls = [[1 + v * k + c for c in range(k)] for v in range(size)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        if sv.solve():
            mo = sv.get_model()
            col = [next(c for c in range(k) if mo[v * k + c] > 0)
                   for v in range(size)]
            sv.delete()
            return phi, imgs, col
        sv.delete()
    return None, None, None


for rots in ((3, 4), (3, 4, 7)):
    U = units(rots, 1)
    r, ints = coords(U)
    phi, imgs, qcol = find_phi(ints, r)
    if phi is None:
        print(f"{rots}: no phi found", flush=True)
        continue
    print(f"\nrotations {rots}: rank {r}, phi to Z_5 x Z_5 found; "
          f"quotient colouring {qcol}  [{time.time()-t0:.0f}s]", flush=True)

    # walk out some real points, carrying their Z^r coordinates
    ZERO = tuple([0] * r)
    basis = IntBasis.covering(U + [Point(K.zero(), K.zero())])
    urows = basis.rows(U)
    dim = basis.dim
    z = np.zeros(2 * dim, dtype=np.int64)
    pos = {z.tobytes(): ZERO}
    rows = [z]
    q = deque([0])
    while q and len(rows) < 4000:
        i = q.popleft()
        base = rows[i]
        bco = pos[base.tobytes()]
        for s, u in enumerate(urows):
            for sign in (1, -1):
                w = base + sign * u
                b = w.tobytes()
                if b in pos:
                    continue
                co = list(bco)
                co[0] = co[0]        # placeholder, real coords below
                nc = tuple(x + sign * y for x, y in zip(bco, ints[s]))
                pos[b] = nc
                rows.append(w)
                q.append(len(rows) - 1)
                if len(rows) >= 4000:
                    break
            if len(rows) >= 4000:
                break
    arr = np.array(rows, dtype=np.int64)
    E = fast_edges_complete(basis, arr)
    col = []
    for w in rows:
        co = pos[w.tobytes()]
        g = tuple(sum(c * phi[t][s] for t, c in enumerate(co)) % A[s]
                  for s in range(len(A)))
        col.append(qcol[g[0] * A[1] + g[1]])
    bad = [(a, b) for a, b in E if col[a] == col[b]]
    print(f"   {len(rows)} real points, {len(E)} exact unit edges, "
          f"{len(set(col))} colours used, {len(bad)} monochromatic"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"   {'VERIFIED' if not bad else 'FAILS'}: the periodic colouring "
          f"is proper on these points", flush=True)
print("\nDONE", flush=True)
