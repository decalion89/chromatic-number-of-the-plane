"""Density against chromatic number, across everything built here.

"Density is the missing ingredient" drove the last stretch, and it came from a
real observation -- Sa and G sit at 4.98 edges per vertex and no operation in
the family moves it.  But there is a result earlier in this same work that
cuts against it: SINGLE-DISTANCE LATTICES ARE BIPARTITE.  The densest
unit-distance graphs known are lattice-like, and lattice-like means
2-chromatic.  So pushing density pushes toward the achromatic end, and the two
goals are in tension rather than aligned.

That is worth measuring rather than asserting, so: for each family built here,
the points, the edges per vertex, and the actual chromatic number computed by
SAT -- smallest k with a proper k-colouring.  If the tension is real the table
will show the densest graphs with the lowest chromatic numbers.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import deque
from hn.degrey import build_Sa, build_Y, build_G
from hn.field import Field
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from hn.fast import IntBasis, fast_edges_complete
from quotient import units
from pysat.solvers import Solver

t0 = time.time()
K = Field((3, 5, 7, 11))
ZERO = Point(K.zero(), K.zero())


def chrom(n, E, hi=6):
    for k in range(2, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        if ok:
            return k
    return f">{hi}"


def from_points(P):
    g = build_graph(P)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    return g.n, E


def walk(U, R, cap):
    basis = IntBasis.covering(U + [ZERO])
    srows = basis.rows(U)
    dim, D = basis.dim, basis.D
    steps = np.concatenate([srows, -srows], axis=0)
    lim = R * R * D * D
    f = lambda V: (basis._field_square(V[:, :dim])
                   + basis._field_square(V[:, dim:])) @ basis.sqrts
    z = np.zeros(2 * dim, dtype=np.int64)
    seen, out, q = {z.tobytes()}, [z], deque([0])
    while q and len(out) < cap:
        i = q.popleft()
        cand = steps + out[i]
        for w in cand[f(cand) <= lim]:
            b = w.tobytes()
            if b not in seen and len(out) < cap:
                seen.add(b)
                out.append(w)
                q.append(len(out) - 1)
    rows = np.array(out, dtype=np.int64)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    return len(rows), E


def tri(R):
    half, rt3 = K.rational(Fr(1, 2)), K.sqrt(3) * K.rational(Fr(1, 2))
    P = []
    n = int(R) + 2
    for a in range(-n, n + 1):
        for b in range(-n, n + 1):
            if a * a + a * b + b * b <= R * R:
                P.append(Point(K.rational(a) + half * K.rational(b),
                               rt3 * K.rational(b)))
    return P


ROWS = []
for name, (n, E) in (
        ("triangular lattice R=12", from_points(tri(12))),
        ("Sa", from_points(build_Sa(F))),
        ("Y", from_points(build_Y(F))),
        # G is skipped: its chi = 5 needs the k=4 UNSAT proof, which is de
        # Grey's own hard result, already reproduced and recorded here.

        ("walk rho_4 e<=2, 7000", walk(units((4,), 2), 2.0, 7000)),
        ("walk rho_3,rho_4, 7000", walk(units((3, 4), 1), 2.0, 7000)),
        ("walk rho_3,4,7, 7000", walk(units((3, 4, 7), 1), 2.0, 7000)),
):
    c = chrom(n, E)
    ROWS.append((name, n, len(E), len(E) / n, c))
    print(f"  {name:28s} {n:6d} pts  {len(E)/n:5.2f} per vertex  chi = {c}"
          f"  [{time.time()-t0:.0f}s]", flush=True)

print("\nsorted by density:", flush=True)
for name, n, e, d, c in sorted(ROWS, key=lambda r: -r[3]):
    print(f"  {d:5.2f} per vertex   chi = {c}   {name} ({n} pts)", flush=True)
print("\nDONE", flush=True)
