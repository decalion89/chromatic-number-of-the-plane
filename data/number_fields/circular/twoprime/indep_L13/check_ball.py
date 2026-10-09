#!/usr/bin/env python3
"""Exact check of the radius-2 ball claim (end of the 'Finite witnesses' subsection).

U = G_25 V over Q(sqrt 11), V = {1, u_n, conj u_n : n = 1, 7, 19}, u_n = (n + i sqrt11)^2/(n^2 + 11),
G_25 = {i^a rho^j : |j| <= 2}, rho = (3+4i)/5.  Points of L = F(i) are stored exactly as integer
4-tuples (X0, X1, Y0, Y1) over a common denominator D:  ((X0 + X1 s) + (Y0 + Y1 s) i)/D, s = sqrt 11.
"""
import sys
from fractions import Fraction as Fr
from math import lcm
from collections import deque
import numpy as np

d = 11
FAILS = []


def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""))
    if not cond:
        FAILS.append(name)


# exact elements of L as 4-tuples of Fractions (x0, x1, y0, y1): (x0 + x1 s) + (y0 + y1 s) i
def lmul(a, b):
    def fmul(p, q):  # (p0 + p1 s)(q0 + q1 s)
        return (p[0] * q[0] + d * p[1] * q[1], p[0] * q[1] + p[1] * q[0])
    ax, ay, bx, by = (a[0], a[1]), (a[2], a[3]), (b[0], b[1]), (b[2], b[3])
    xx, yy, xy, yx = fmul(ax, bx), fmul(ay, by), fmul(ax, by), fmul(ay, bx)
    return (xx[0] - yy[0], xx[1] - yy[1], xy[0] + yx[0], xy[1] + yx[1])


def lconj(a):
    return (a[0], a[1], -a[2], -a[3])


def lnorm(a):  # a * conj(a) as an element of F
    p = lmul(a, lconj(a))
    assert p[2] == 0 and p[3] == 0
    return (p[0], p[1])


ONE = (Fr(1), Fr(0), Fr(0), Fr(0))
I = (Fr(0), Fr(0), Fr(1), Fr(0))
RHO = (Fr(3, 5), Fr(0), Fr(4, 5), Fr(0))
RHOI = lconj(RHO)


def lpow(a, n):
    r = ONE
    for _ in range(n):
        r = lmul(r, a)
    return r


G25 = []
for a in range(4):
    for j in range(-2, 3):
        G25.append(lmul(lpow(I, a), lpow(RHO, j) if j >= 0 else lpow(RHOI, -j)))
check("|G_25| = 20 distinct rational rotations of norm 1",
      len(set(G25)) == 20 and all(lnorm(g) == (1, 0) for g in G25))

V = [ONE]
for n in (1, 7, 19):
    un = lmul(lmul((Fr(n), Fr(0), Fr(0), Fr(1)), (Fr(n), Fr(0), Fr(0), Fr(1))),
              (Fr(1, n * n + d), Fr(0), Fr(0), Fr(0)))
    V += [un, lconj(un)]
check("V: 7 distinct unit vectors", len(set(V)) == 7 and all(lnorm(v) == (1, 0) for v in V))
U = sorted({lmul(g, v) for g in G25 for v in V})
check("U = G_25 V has 140 elements, symmetric, all of norm 1",
      len(U) == 140 and all(tuple(-x for x in u) in set(U) for u in U) and all(lnorm(u) == (1, 0) for u in U))

D = 1
for u in U:
    for x in u:
        D = lcm(D, x.denominator)
print("common denominator D =", D)
Uint = [tuple(int(x * D) for x in u) for u in U]


def iadd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


ZERO = (0, 0, 0, 0)
layer = {ZERO: 0}
for u in Uint:
    layer.setdefault(u, 1)
for u in Uint:
    for v in Uint:
        layer.setdefault(iadd(u, v), 2)
pts = list(layer)
idx = {p: k for k, p in enumerate(pts)}
n = len(pts)
check("ball of radius 2 has 9941 vertices", n == 9941, "n = %d, layers %s" % (
    n, [sum(1 for p in pts if layer[p] == t) for t in (0, 1, 2)]))

# Cayley edges inside the ball
cay = set()
for p in pts:
    for u in Uint:
        q = iadd(p, u)
        if q in idx:
            a, b = idx[p], idx[q]
            cay.add((min(a, b), max(a, b)))
check("the induced Cayley graph has 19600 edges", len(cay) == 19600, "%d" % len(cay))
same_layer = sum(1 for (a, b) in cay if layer[pts[a]] == layer[pts[b]])
check("no Cayley edge joins two vertices of the same layer", same_layer == 0)

adjc = [[] for _ in range(n)]
for a, b in cay:
    adjc[a].append(b)
    adjc[b].append(a)


def bipartite(adj):
    col = [-1] * n
    for s0 in range(n):
        if col[s0] >= 0:
            continue
        col[s0] = 0
        dq = deque([s0])
        while dq:
            x = dq.popleft()
            for y in adj[x]:
                if col[y] < 0:
                    col[y] = 1 - col[x]
                    dq.append(y)
                elif col[y] == col[x]:
                    return False
    return True


check("the Cayley ball is bipartite", bipartite(adjc))

# all unit-distance pairs, exactly:  (dX0 + dX1 s)^2 + (dY0 + dY1 s)^2 = D^2
A = np.array(pts, dtype=np.int64)
assert np.abs(A).max() < 10 ** 6
D2 = D * D
unit = []
for k in range(n - 1):
    dd = A[k + 1:] - A[k]
    r = dd[:, 0] * dd[:, 0] + d * dd[:, 1] * dd[:, 1] + dd[:, 2] * dd[:, 2] + d * dd[:, 3] * dd[:, 3]
    s_ = dd[:, 0] * dd[:, 1] + dd[:, 2] * dd[:, 3]
    for m in np.nonzero((r == D2) & (s_ == 0))[0]:
        unit.append((k, k + 1 + int(m)))
unit = set(unit)
check("every Cayley edge is a unit-distance pair", cay <= unit)
extra = unit - cay
check("the induced unit-distance graph has 216 more edges", len(extra) == 216, "%d extra" % len(extra))
from collections import Counter
print("layers of the extra edges:", Counter(tuple(sorted((layer[pts[a]], layer[pts[b]]))) for a, b in extra))
# exact re-verification of the extra edges with Fractions
ok_exact = True
for a, b in extra:
    diff = tuple(Fr(pts[b][t] - pts[a][t], D) for t in range(4))
    if lnorm(diff) != (1, 0) or tuple(int(x * D) for x in diff) in set(Uint):
        ok_exact = False
check("each extra edge: difference has norm exactly 1 and is not in U", ok_exact)

adj = [set() for _ in range(n)]
for a, b in unit:
    adj[a].add(b)
    adj[b].add(a)
tri = sum(1 for a, b in unit if adj[a] & adj[b])
check("the unit-distance graph has no triangle", tri == 0)

# a 5-cycle (it must use an extra edge, as the Cayley part is bipartite)
# search properly for a 5-cycle a-x-z-y-b-a
cyc5 = None
for a, b in sorted(extra):
    for x in adj[a] - {b}:
        for y in adj[b] - {a, x}:
            common = (adj[x] & adj[y]) - {a, b}
            if common:
                z = min(common)
                cyc5 = (a, x, z, y, b)
                break
        if cyc5:
            break
    if cyc5:
        break
ok5 = cyc5 is not None and len(set(cyc5)) == 5 and all(cyc5[(t + 1) % 5] in adj[cyc5[t]] for t in range(5))
check("the unit-distance graph contains a 5-cycle", ok5, str([layer[pts[v]] for v in cyc5]) if cyc5 else "")

# homomorphism to K_{5/2} = C_5 (colours Z/5, adjacent iff difference is 2 or 3 mod 5)
# reduction: c(0) = 0, c(U) in {2,3}; a layer-2 vertex with no extra edge can take colour 0.
OK = lambda x, y: (y - x) % 5 in (2, 3)
core = sorted({v for e in extra for v in e} | {idx[ZERO]} | {v for e in extra for v in e for w in adj[v] if layer[pts[w]] == 1}
              | {w for e in extra for v in e for w in adj[v] if layer[pts[w]] == 1})
core_set = set(core)
dom = {}
for v in core:
    L_ = layer[pts[v]]
    dom[v] = [0] if L_ == 0 else ([2, 3] if L_ == 1 else [0, 1, 4])
col = {}


def cands(v):
    return [c for c in dom[v] if all(OK(col[w], c) for w in adj[v] if w in col)]


nodes = [0]


def solve():
    free = [v for v in core if v not in col]
    if not free:
        return True
    nodes[0] += 1
    if nodes[0] > 2_000_000:
        raise RuntimeError("search limit")
    v = min(free, key=lambda x: (len(cands(x)), -len(adj[x] & core_set)))
    for c in cands(v):
        col[v] = c
        if solve():
            return True
        del col[v]
    return False


sys.setrecursionlimit(100000)
found = solve()
print('search nodes:', nodes[0])
full = None
if found:
    full = [None] * n
    for v in range(n):
        L_ = layer[pts[v]]
        if v in col:
            full[v] = col[v]
        elif L_ == 1:
            full[v] = 2
        else:
            full[v] = 0
ok_hom = full is not None and all(OK(full[a], full[b]) for a, b in unit)
check("the unit-distance graph has a homomorphism to K_{5/2} (verified on all %d edges)" % len(unit), ok_hom,
      "core size %d" % len(core))

print()
print("FAILED:", FAILS if FAILS else "none")
sys.exit(1 if FAILS else 0)
