"""mkcnf_ball.py NAME -- write the 4-colouring CNF (triangle 0,1,zeta6 fixed) for a ball test graph.
NAME: BDk = B_k(D) u conj(B_k(D)); BUk = B_k(U); the 4-core is taken first (vertices of degree < 4 removed
iteratively, which does not change 4-colourability)."""
import sys, time
import numpy as np
from flat import *
from satutil import write_dimacs

D, CD, U = directions()
name = sys.argv[1]
k = int(name[2:])
if name.startswith("BD"):
    B = sums_ball(D, k)
    P = unique_rows(np.concatenate([B, conj_rows(B)]))
elif name.startswith("BU"):
    P = sums_ball(U, k)
n0 = len(P)
E, J = build_edges(P, U)
# 4-core
alive = np.ones(len(P), bool)
while True:
    deg = np.bincount(E[alive[E[:, 0]] & alive[E[:, 1]]].ravel(), minlength=len(P))
    low = alive & (deg < 4)
    if not low.any():
        break
    alive &= ~low
P = P[alive]
E, J = build_edges(P, U)
S = PointSet(P)
tri = tuple(int(S.index([x])[0]) for x in ([0] * 12, [7] + [0] * 11, D[14]))
assert min(tri) >= 0
n = len(P)
cls = cnf_clauses(n, E.tolist(), tri)
write_dimacs(f"{name}.cnf", 4 * n, cls, [f"{name}: 4-core of {n0} points -> {n} vertices, {len(E)} edges; triangle {tri} fixed to colours 1,2,3"])
np.save(f"pts_{name}.npy", P)
print(f"{name}: {n0} points, 4-core {n} vertices, {len(E)} edges, {len(cls)} clauses")
