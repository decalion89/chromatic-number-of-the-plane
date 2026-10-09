"""Referee's own exact enumeration of the components of S^{r0}(K,M) modulo N = 5^K 13^M by convex-polygon clipping
(no code taken from the authors' programs).

For each cell [a+r0, a+1-r0] x [b+r0, b+1-r0] (0 <= a, b < N) the square is clipped successively by the strips
n + r0 <= f <= n + 1 - r0 of every other functional f (Re and Im of conj(c) rho^j sigma^l); all integers n whose strip
meets the current polygon are followed.  Nonempty leaves give the index vectors.  Exact Fractions throughout.

usage: python3 clip_window.py K M r0 certfile
"""
from fractions import Fraction as Fr
import sys
import time


def cdiv(x, y):
    n = y[0] * y[0] + y[1] * y[1]
    return (Fr(x[0] * y[0] + x[1] * y[1], n), Fr(x[1] * y[0] - x[0] * y[1], n))


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z, e = cdiv((1, 0), z), -e
    for _ in range(e):
        w = cmul(w, z)
    return w


def fl(q):
    return q.numerator // q.denominator


K, M = int(sys.argv[1]), int(sys.argv[2])
r0 = Fr(sys.argv[3])
certfile = sys.argv[4]
N = 5 ** K * 13 ** M
RHO, SIG = cdiv((2, 1), (2, -1)), cdiv((3, 2), (3, -2))
KEYS = [(j, l, e) for j in range(-K, K + 1) for l in range(-M, M + 1) for e in ("Re", "Im")]
FUN = {}
for (j, l, e) in KEYS:
    g = cmul(cpow(RHO, j), cpow(SIG, l))
    FUN[(j, l, e)] = (g[0], g[1]) if e == "Re" else (g[1], -g[0])
# conj(c) g with c = x + iy: Re = g0 x + g1 y, Im = g1 x - g0 y
assert FUN[(0, 0, "Re")] == (1, 0) and FUN[(0, 0, "Im")] == (0, -1)
OTHER = [k for k in KEYS if k[:2] != (0, 0)]
# integer scaling of each functional (all coefficients have denominator dividing N)
LO, HI = r0, 1 - r0


def clip(poly, a, b, c):
    """keep a x + b y <= c (closed); poly is a list of vertices of a convex polygon (possibly degenerate)"""
    out = []
    n = len(poly)
    for i in range(n):
        P, Q = poly[i], poly[(i + 1) % n]
        vp = a * P[0] + b * P[1] - c
        vq = a * Q[0] + b * Q[1] - c
        if vp <= 0:
            out.append(P)
        if (vp < 0 < vq) or (vq < 0 < vp):
            t = vp / (vp - vq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    # remove consecutive duplicates
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    if len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


found = {}
nodes = 0


def rec(depth, poly, idx):
    global nodes
    if depth == len(OTHER):
        iv = tuple(idx[k] for k in KEYS)
        found[iv] = found.get(iv, 0) + 1
        return
    k = OTHER[depth]
    fa, fb = FUN[k]
    vals = [fa * p[0] + fb * p[1] for p in poly]
    mn, mx = min(vals), max(vals)
    for n in range(fl(mn - HI), fl(mx - LO) + 1):
        if n + HI < mn or n + LO > mx:
            continue
        p2 = clip(poly, fa, fb, n + HI)
        if not p2:
            continue
        p2 = clip(p2, -fa, -fb, -(n + LO))
        if not p2:
            continue
        nodes += 1
        idx[k] = n
        rec(depth + 1, p2, idx)
        del idx[k]


t0 = time.time()
for a in range(N):
    for b in range(N):
        sq = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
        rec(0, sq, {(0, 0, "Re"): a, (0, 0, "Im"): -b - 1})
print(f"(K,M)=({K},{M}) N={N} r0={r0}: {len(found)} index vectors with nonempty polygons, {nodes} nonempty branch "
      f"nodes, {time.time() - t0:.0f}s")

# compare with the certificate's list (same functional names)
lines = open(certfile).read().splitlines()
fl_line = [ln for ln in lines if ln.startswith("# functionals")][0]
ckeys = [tuple(x.strip("()").split(",")) for x in fl_line.split(":", 1)[1].split()]
ckeys = [(int(j), int(l), e) for (j, l, e) in ckeys]
cert = {}
for ln in lines:
    if ln.startswith("component"):
        head, lab = ln.split("|", 1)
        v = [int(x) for x in head.split("indices")[1].split()]
        d = dict(zip(ckeys, v))
        cert[tuple(d[k] for k in KEYS)] = lab.strip().split()[0]
print("same set of index vectors as the certificate:", set(found) == set(cert),
      "| labels:", {lab: sum(1 for v in cert.values() if v == lab) for lab in set(cert.values())})
