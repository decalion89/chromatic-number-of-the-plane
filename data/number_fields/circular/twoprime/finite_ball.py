"""The ball of radius 2 in the Cayley graph of the 140 unit vectors of Theorem C over Q(sqrt11) (paper Section 10,
"Finite witnesses"; note section 6.8), with exact arithmetic.

U = G_25 V, V = {1, u_n, conj u_n : n = 1, 7, 19}, u_n = (n + i sqrt11)^2/(n^2 + 11), G_25 = {i^a rho^j : |j| <= 2},
rho = (3 + 4i)/5. The ball B = {0} u U u (U + U) has 9 941 points. The program checks:
(1) the Cayley edges (pairs whose difference lies in U) number 19 600, and they form a bipartite graph;
(2) the unit-distance graph induced on B has 216 further edges, no triangle, and a 5-cycle (printed);
(3) it maps to K_{5/2} = C_5 (a colouring found by SAT and rechecked), so its circular chromatic number is 5/2.
usage: python3 finite_ball.py"""
from collections import deque
from fractions import Fraction as Q
from math import lcm
import numpy as np
from pysat.solvers import Cadical153

d = 11
# elements of L = F(i), F = Q(sqrt d), as (a, b, c, e): x = a + b r, y = c + e r, z = x + i y, r = sqrt d


def lmul(z, w):
    def fmul(p, q):
        return (p[0] * q[0] + d * p[1] * q[1], p[0] * q[1] + p[1] * q[0])
    x1, y1, x2, y2 = (z[0], z[1]), (z[2], z[3]), (w[0], w[1]), (w[2], w[3])
    a = fmul(x1, x2); b = fmul(y1, y2); c = fmul(x1, y2); e = fmul(y1, x2)
    return (a[0] - b[0], a[1] - b[1], c[0] + e[0], c[1] + e[1])


def u(n):
    D = n * n + d
    return (Q(n * n - d, D), Q(0), Q(0), Q(2 * n, D))


ONE, I = (Q(1), Q(0), Q(0), Q(0)), (Q(0), Q(0), Q(1), Q(0))
V = [ONE] + [w for n in (1, 7, 19) for w in (u(n), (u(n)[0], u(n)[1], -u(n)[2], -u(n)[3]))]
rho, rhoi = (Q(3, 5), Q(0), Q(4, 5), Q(0)), (Q(3, 5), Q(0), Q(-4, 5), Q(0))
G = []
for a in range(4):
    ia = ONE
    for _ in range(a):
        ia = lmul(ia, I)
    for j in range(-2, 3):
        g = ia
        for _ in range(abs(j)):
            g = lmul(g, rho if j > 0 else rhoi)
        G.append(g)
U = sorted({lmul(g, v) for g in G for v in V})
assert len(U) == 140
for z in U:   # x^2 + y^2 = 1 in F, exactly
    assert (z[0] ** 2 + d * z[1] ** 2 + z[2] ** 2 + d * z[3] ** 2, 2 * z[0] * z[1] + 2 * z[2] * z[3]) == (1, 0)
Uset = set(U)
zero = (Q(0),) * 4
ball = {zero} | Uset | {tuple(p[k] + w[k] for k in range(4)) for p in U for w in U}
pts = sorted(ball)
n = len(pts)
M = lcm(*(t.denominator for p in pts for t in p))
P = [tuple(int(t * M) for t in p) for p in pts]
A = np.array(P, dtype=np.int64)
assert all(int(A[i, k]) == P[i][k] for i in range(n) for k in range(4))
edges = []
for i in range(n):
    dd = A[i + 1:] - A[i]
    rr = dd[:, 0] * dd[:, 0] + d * dd[:, 1] * dd[:, 1] + dd[:, 2] * dd[:, 2] + d * dd[:, 3] * dd[:, 3]
    ss = dd[:, 0] * dd[:, 1] + dd[:, 2] * dd[:, 3]
    for j in np.nonzero((rr == M * M) & (ss == 0))[0]:
        edges.append((i, i + 1 + int(j)))
for i, j in edges:   # exact recheck with Python integers
    a = [P[j][k] - P[i][k] for k in range(4)]
    assert a[0] ** 2 + d * a[1] ** 2 + a[2] ** 2 + d * a[3] ** 2 == M * M and a[0] * a[1] + a[2] * a[3] == 0
cay = [(i, j) for i, j in edges if tuple(pts[j][k] - pts[i][k] for k in range(4)) in Uset]
print(f"ball of radius 2: {n} points; unit-distance pairs {len(edges)}; Cayley edges {len(cay)}; "
      f"others {len(edges) - len(cay)}")


def adjacency(E):
    adj = [[] for _ in range(n)]
    for i, j in E:
        adj[i].append(j); adj[j].append(i)
    return adj


def two_colour(E):
    adj, col = adjacency(E), [-1] * n
    for s in range(n):
        if col[s] < 0:
            col[s] = 0; q = deque([s])
            while q:
                x = q.popleft()
                for y in adj[x]:
                    if col[y] < 0:
                        col[y] = 1 - col[x]; q.append(y)
                    elif col[y] == col[x]:
                        return None
    return col


col = two_colour(cay)
assert col is not None and all(col[i] != col[j] for i, j in cay)
print("Cayley ball bipartite: True (2-colouring rechecked)")
assert two_colour(edges) is None
adj = adjacency(edges)
nb = [set(a) for a in adj]
assert sum(len(nb[i] & nb[j]) for i, j in edges) == 0
print("induced unit-distance graph: not bipartite, no triangle")
# a 5-cycle x - a - b - y - c - x (five distinct vertices)
cyc = None
for x in range(n):
    if cyc:
        break
    for a in nb[x]:
        if cyc:
            break
        for c in nb[x]:
            if c <= a:
                continue
            # need b adjacent to a, y adjacent to c, and b ~ y, all distinct
            for b in nb[a]:
                if b in (x, c):
                    continue
                ys = (nb[b] & nb[c]) - {x, a}
                if ys:
                    cyc = [x, a, b, min(ys), c]
                    break
            if cyc:
                break
assert cyc and len(set(cyc)) == 5 and all(cyc[(k + 1) % 5] in nb[cyc[k]] for k in range(5))
print("5-cycle:", cyc)
var = lambda v, c: 5 * v + c + 1
s = Cadical153()
for v in range(n):
    s.add_clause([var(v, c) for c in range(5)])
    for a in range(5):
        for b in range(a + 1, 5):
            s.add_clause([-var(v, a), -var(v, b)])
for i, j in edges:
    for a in range(5):
        for b in range(5):
            if (b - a) % 5 not in (2, 3):
                s.add_clause([-var(i, a), -var(j, b)])
assert s.solve()
m = set(x for x in s.get_model() if x > 0)
c5 = [next(c for c in range(5) if var(v, c) in m) for v in range(n)]
assert all((c5[j] - c5[i]) % 5 in (2, 3) for i, j in edges)
print("homomorphism to K_{5/2} = C_5: found and rechecked; circular chromatic number 5/2")
