"""End-to-end exact test of Lemma 14 (lem:window7) at level k = 6 (and of Lemma 12 for theta > 2/7).

S^theta(k,1) is computed exactly as a union of convex rational polygons modulo P_k = 13*5^k:
level 0 (rotations i^a sigma^l) by clipping the 13^2 cells, then level k from level k-1 by the 25 translates by
P_{k-1}(u+vi) and clipping with the 12 new functionals Re/Im of conj(c) rho^{+-k} sigma^l.

At the last level every component (a polygon) is tested against the conclusion of Lemma 14, with the text's own
definitions: for some eps in E' (13 classes mod Z[i]) and m in Z[i], every vertex v gives x = v - N eps - N m with
  type c : conj(x) rho^j in B_s                      (|j| <= k)   [x in P_k]
  type q : eta_j(eps)/6 + conj(x) rho^j in B_s        (|j| <= k)   [x in Y_k(eps)]
  type 7 : zeta_j(eps) + conj(x) rho^j in B_s         (|j| <= k)   [x in X_k(eps)], zeta_j via rho~ = 2+5i mod 7
and |x| <= r' = (1/4+3/14) sqrt 2 (checked as |x|^2 <= 2 (13/28)^2).  Polygons are convex and the sets are convex,
so vertices suffice.

usage: python3 lift_lemma14.py theta K
"""
from fractions import Fraction as Fr
import sys
import time

theta = Fr(sys.argv[1])
KMAX = int(sys.argv[2])
s = Fr(1, 2) - theta
LO, HI = theta, 1 - theta


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


def conj(x):
    return (x[0], -x[1])


def fl(q):
    return q.numerator // q.denominator


RHO, SIG = cdiv((2, 1), (2, -1)), cdiv((3, 2), (3, -2))


def funcs(j, l):
    g = cmul(cpow(RHO, j), cpow(SIG, l))
    return [(g[0], g[1]), (g[1], -g[0])]     # Re and Im of conj(c) g, c = x + iy


def clip(poly, a, b, c):
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
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    if len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


def cut(pieces, fs):
    """clip every piece by the strips n+theta <= f <= n+1-theta of every functional in fs (all n)"""
    for (fa, fb) in fs:
        new = []
        for poly in pieces:
            vals = [fa * p[0] + fb * p[1] for p in poly]
            mn, mx = min(vals), max(vals)
            for n in range(fl(mn - HI), fl(mx - LO) + 1):
                p2 = clip(poly, fa, fb, n + HI)
                if p2:
                    p2 = clip(p2, -fa, -fb, -(n + LO))
                if p2:
                    new.append(p2)
        pieces = new
    return pieces


t0 = time.time()
# level 0: rotations i^a sigma^l, |l| <= 1, period 13
P0 = 13
F0 = funcs(0, 1) + funcs(0, -1)
comps = []
for a in range(P0):
    for b in range(P0):
        sq = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
        comps += cut([sq], F0)
print(f"theta={theta}: level 0, period 13: {len(comps)} components ({time.time() - t0:.1f}s)")
per = P0
for k in range(1, KMAX + 1):
    Fk = []
    for l in (-1, 0, 1):
        Fk += funcs(k, l) + funcs(-k, l)
    new = []
    for poly in comps:
        for u in range(5):
            for v in range(5):
                tr = [(p[0] + per * u, p[1] + per * v) for p in poly]
                new += cut([tr], Fk)
    comps = new
    per *= 5
    print(f"  level {k}, period {per}: {len(comps)} components ({time.time() - t0:.1f}s)")

# ------------------------------------------------------------------ test the conclusion of Lemma 14 at level KMAX
k = KMAX
N = 5 ** k
H = (Fr(1, 2), Fr(1, 2))
EC = [H]
EQ = [(Fr(a, 3), Fr(b, 3)) for a in (1, 2) for b in (1, 2)]
A7 = [(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)]
E7 = [(Fr(a, 7), Fr(b, 7)) for (a, b) in A7]
RHOJ = {j: cpow(RHO, j) for j in range(-k, k + 1)}


def inBs(w):
    return -s <= w[0] <= s and -s <= w[1] <= s


def gpart(w):
    return (w[0] - H[0] - fl(w[0]), w[1] - H[1] - fl(w[1]))


def eta6(eps, j):
    # eta_j(eps)/6 from the definition conj(N eps) rho^j in h + eta/6 + Z[i]
    g = gpart(cmul(conj((N * eps[0], N * eps[1])), RHOJ[j]))
    assert abs(g[0]) == Fr(1, 6) and abs(g[1]) == Fr(1, 6)
    return g


def m7(x, y):
    return ((x[0] * y[0] - x[1] * y[1]) % 7, (x[0] * y[1] + x[1] * y[0]) % 7)


RT = [(1, 0)]
for _ in range(7):
    RT.append(m7(RT[-1], (2, 5)))


def zeta(eps, j):
    # conj(eps) rho~^j in h + zeta + Z[i], rho~ = 2+5i, rho~^j for j<0 via rho~^8 = 1 mod 7
    e = (int(7 * eps[0]), int(7 * eps[1]))
    w = m7((e[0], -e[1]), RT[j % 8])
    g = gpart((Fr(w[0], 7), Fr(w[1], 7)))
    assert all(abs(t) in (Fr(1, 14), Fr(3, 14)) for t in g)
    return g


R2 = 2 * Fr(13, 28) ** 2            # r'^2


def fits(poly, eps, kind):
    v0 = poly[0]
    m = (round((v0[0] - N * eps[0]) / N), round((v0[1] - N * eps[1]) / N))
    off = {j: (eta6(eps, j) if kind == "q" else zeta(eps, j) if kind == "7" else (Fr(0), Fr(0)))
           for j in range(-k, k + 1)}
    worst = Fr(0)
    for v in poly:
        x = (v[0] - N * (eps[0] + m[0]), v[1] - N * (eps[1] + m[1]))
        for j in range(-k, k + 1):
            w = cmul(conj(x), RHOJ[j])
            if not inBs((off[j][0] + w[0], off[j][1] + w[1])):
                return None
        worst = max(worst, x[0] ** 2 + x[1] ** 2)
    return worst


count = {"c": 0, "q": 0, "7": 0, "none": 0}
worst = {"c": Fr(0), "q": Fr(0), "7": Fr(0)}
for poly in comps:
    hit = None
    for kind, lst in (("c", EC), ("q", EQ), ("7", E7)):
        for eps in lst:
            w = fits(poly, eps, kind)
            if w is not None:
                hit = (kind, w)
                break
        if hit:
            break
    if hit is None:
        count["none"] += 1
    else:
        count[hit[0]] += 1
        worst[hit[0]] = max(worst[hit[0]], hit[1])
print(f"level {k} (N = 5^{k} = {N % 7} mod 7): components by the form of Lemma 14 / Lemma 12: {count}")
print("  largest |x|^2 by type:", {t: float(w) for t, w in worst.items()},
      "; all <= r'^2 =", float(R2), ":", all(w <= R2 for w in worst.values()))
