"""The same confined ball, in integers -- the exact version was too slow.

Sums of unit vectors do not grow denominators: w has denominator two, rho_4
has eight, so every product w^a rho_4^e has denominator dividing sixteen and
every SUM of them does too.  The whole construction therefore lives in one
fixed lattice, and can be done in int64 with numpy instead of in Fractions
over a sixteen-dimensional field, which is where the first attempt stalled --
twelve minutes and not one ball finished.

Everything else is unchanged.  Unit vectors are products w^a rho_D^e, so a few
rotations give dozens of directions.  Points are kept only inside a small disc,
so that a point's neighbours are kept too and the degree survives rather than
bleeding away at a boundary -- which is what made the outward-growing Cayley
ball come in at 4.93 per vertex against the sumset's 6.59.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import deque
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from hn.graph import build_graph
from hn.fast import IntBasis, fast_edges_complete
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


def confined(U, R, cap=20000):
    """BFS in int64: rows are numerators over the shared denominator."""
    basis = IntBasis.covering(U + [Point(K.zero(), K.zero())])
    rows = basis.rows(U)
    D = basis.D
    steps = np.concatenate([rows, -rows], axis=0)
    dim = basis.dim
    lim = R * R * D * D

    # the norms of all candidates in one call, not one call each -- the
    # per-candidate version was the whole cost
    def norms(V):
        sq = (basis._field_square(V[:, :dim])
              + basis._field_square(V[:, dim:]))
        return sq @ basis.sqrts

    z = np.zeros(2 * dim, dtype=np.int64)
    seen = {z.tobytes()}
    out = [z]
    q = deque([z])
    while q and len(out) < cap:
        p = q.popleft()
        cand = steps + p
        good = cand[norms(cand) <= lim]
        for w in good:
            b = w.tobytes()
            if b in seen:
                continue
            seen.add(b)
            out.append(w)
            q.append(w)
            if len(out) >= cap:
                break
    return basis, np.array(out, dtype=np.int64)


def solve(basis, rows):
    """Edges from the int64 rows directly.

    build_graph on exact Points is O(n^2) in Fraction arithmetic and stalls
    past a few thousand; fast_edges_complete buckets into cells of side one,
    nominates candidates by float distance and confirms each by exact integer
    field arithmetic, so the edges are still exactly the unit-distance pairs.
    """
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    n = len(rows)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return n, len(E), ok


for rots, emax in (((4,), 1), ((3, 4), 1), ((4,), 2), ((3, 4, 7), 1)):
    U = units(rots, emax)
    print(f"\nrotations {rots}, exponents <= {emax}: {len(U)} unit vectors"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for R in (1.5, 2.0, 2.5):
        basis, rows = confined(U, R)
        if len(rows) < 60:
            print(f"   R={R}: only {len(rows)} points"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            continue
        print(f"   R={R}: BFS gave {len(rows)} points, building the graph"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        n, e, ok = solve(basis, rows)
        print(f"   R={R}: {n} pts, {e} edges, {e/n:.2f} per vertex -> "
              f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if not ok:
            print("*** chi >= 6 ***", flush=True)
            sys.exit(0)
print("\nDONE", flush=True)
