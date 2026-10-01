"""Periodic 3-colourings of the unit-distance graph over Q(sqrt d) with directions U_D:
is Cay(M/mM, U) 3-colourable, where M is the Z-span of U_D (a lattice in Z^4)?
A yes for some m gives a 3-colouring of every graph built from U_D (the growth cannot succeed).
usage: periodicq.py d D m1 m2 ..."""
import sys, itertools, time
import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from pysat.solvers import Solver
from units_fast import units_fast
d, D = int(sys.argv[1]), int(sys.argv[2])
ms = [int(x) for x in sys.argv[3:]]
U = [tuple(int(x) for x in u) for u in np.asarray(units_fast(d, D)).reshape(-1, 4)]
H = hermite_normal_form(Matrix(U).T)          # columns span the same lattice as the u's
B = Matrix.hstack(*[H[:, j] for j in range(H.shape[1]) if any(H[:, j])])
r = B.shape[1]
Binv = B.inv()
C = [tuple(int(x) for x in Binv * Matrix(u)) for u in U]
assert all(B * Matrix(c) == Matrix(u) for c, u in zip(C, U))
print(f"d={d} D={D}: {len(U)} directions, rank {r}, index of M in Z^4: {abs(B.det()) if r == 4 else 'n/a'}", flush=True)
for m in ms:
    t = time.time()
    G = sorted(set(tuple(x % m for x in c) for c in C))
    if tuple([0] * r) in G:
        print(f"  m={m}: a direction is 0 mod m: no periodic colouring", flush=True)
        continue
    pts = list(itertools.product(range(m), repeat=r))
    idx = {q: i for i, q in enumerate(pts)}
    var = lambda i, c: 3 * i + c + 1
    s = Solver(name='cd19')
    for i in range(len(pts)):
        s.add_clause([var(i, c) for c in range(3)])
    s.add_clause([var(0, 0)])
    for i, q in enumerate(pts):
        for g in G:
            j = idx[tuple((a + b) % m for a, b in zip(q, g))]
            if i < j:
                for c in range(3):
                    s.add_clause([-var(i, c), -var(j, c)])
    ok = s.solve()
    print(f"  m={m}: {len(pts)} classes, {len(G)} generators mod m: {'3-COLOURABLE' if ok else 'not 3-colourable'}  [{time.time() - t:.1f}s]", flush=True)
