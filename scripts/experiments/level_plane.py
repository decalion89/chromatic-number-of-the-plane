"""Cay((Z/q^k)^2, U_k): the unit-distance graph of the q-adic plane modulo q^k.

For a real field with an unramified degree-1 place over a prime q = 3 (mod 4), reduction of the
coordinates modulo q^k maps its plane homomorphically to this graph, so its chromatic number bounds
the field's. Writes a DIMACS formula for proper c-colourings, with vertex 0 and a unit triangle
through it pinned to colours 0, 1, 2 when such a triangle exists.

usage: level_plane.py q k c   (writes level_q<q>_k<k>_c<c>.cnf in the current directory)
At q = 11: k = 1 is 5- but not 4-colourable; k = 2 (14 641 vertices) is not 4-colourable either.
"""
import sys, itertools
q, k, c = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
m = q ** k
U = [(a, b) for a in range(m) for b in range(m) if (a * a + b * b - 1) % m == 0]
print(f"q={q} k={k} modulus {m}: |U| = {len(U)} (expected {(q + 1) * q ** (k - 1)} for q = 3 mod 4)", file=sys.stderr)
Uset = set(U)
n = m * m
idx = lambda a, b: (a % m) * m + (b % m)
edges = set()
for a in range(m):
    for b in range(m):
        v = idx(a, b)
        for (u1, u2) in U:
            w = idx(a + u1, b + u2)
            if v < w:
                edges.add((v, w))
print(f"vertices {n}, edges {len(edges)}, degree {len(U)}", file=sys.stderr)
# a unit triangle 0, u, u' with u - u' in U
tri = next(((u, w) for u in U for w in U if u != w and ((u[0] - w[0]) % m, (u[1] - w[1]) % m) in Uset), None)
X = lambda v, col: 1 + v * c + col
cls = [[X(v, col) for col in range(c)] for v in range(n)]
for v, w in edges:
    for col in range(c):
        cls.append([-X(v, col), -X(w, col)])
pins = []
if tri:
    (u, w) = tri
    pins = [(0, 0), (idx(*u), 1), (idx(*w), 2)]
    for v, col in pins:
        cls.append([X(v, col)])
print(f"pinned triangle: {pins}", file=sys.stderr)
with open(f"level_q{q}_k{k}_c{c}.cnf", "w") as f:
    f.write(f"p cnf {n * c} {len(cls)}\n")
    for cl in cls:
        f.write(" ".join(map(str, cl)) + " 0\n")
print(f"clauses {len(cls)}", file=sys.stderr)
