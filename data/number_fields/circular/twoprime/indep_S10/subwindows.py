"""Context claims after Proposition 8: other windows.  Components of {c : Re(conj(c) g) in [theta,1-theta]+Z for g in
i^a W} modulo P (by clipping every cell), and the exact largest least margin of each component (3-variable LP solved
by enumerating triples of tight constraints, exact Fractions).

usage: python3 subwindows.py NAME theta
  NAME = cross  : W = {1, rho^+-1, sigma^+-1}, P = 65
  NAME = p13_17 : W = {sigma^l tau^m : |l|,|m| <= 1}, tau = (4+i)/(4-i), P = 221
  NAME = p5_17  : W = {rho^j tau^m : |j|,|m| <= 1}, P = 85
"""
from fractions import Fraction as Fr
from itertools import combinations
import sys


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


RHO, SIG, TAU = cdiv((2, 1), (2, -1)), cdiv((3, 2), (3, -2)), cdiv((4, 1), (4, -1))
name, theta = sys.argv[1], Fr(sys.argv[2])
if name == "cross":
    W, P = [cpow(RHO, 1), cpow(RHO, -1), cpow(SIG, 1), cpow(SIG, -1)], 65
elif name == "p13_17":
    W, P = [cmul(cpow(SIG, l), cpow(TAU, m)) for l in (-1, 0, 1) for m in (-1, 0, 1) if (l, m) != (0, 0)], 221
elif name == "w21":
    W, P = [cmul(cpow(RHO, l), cpow(SIG, m)) for l in (-2, -1, 0, 1, 2) for m in (-1, 0, 1) if (l, m) != (0, 0)], 325
elif name == "w11":
    W, P = [cmul(cpow(RHO, l), cpow(SIG, m)) for l in (-1, 0, 1) for m in (-1, 0, 1) if (l, m) != (0, 0)], 65
else:
    W, P = [cmul(cpow(RHO, l), cpow(TAU, m)) for l in (-1, 0, 1) for m in (-1, 0, 1) if (l, m) != (0, 0)], 85
for g in W:
    assert (P * g[0]).denominator == 1 and (P * g[1]).denominator == 1
FUN = []
for g in W:
    FUN += [(g[0], g[1]), (g[1], -g[0])]
ALLF = [(Fr(1), Fr(0)), (Fr(0), Fr(-1))] + FUN
LO, HI = theta, 1 - theta


def clip(poly, a, b, c):
    out = []
    n = len(poly)
    for i in range(n):
        Pp, Q = poly[i], poly[(i + 1) % n]
        vp = a * Pp[0] + b * Pp[1] - c
        vq = a * Q[0] + b * Q[1] - c
        if vp <= 0:
            out.append(Pp)
        if (vp < 0 < vq) or (vq < 0 < vp):
            t = vp / (vp - vq)
            out.append((Pp[0] + t * (Q[0] - Pp[0]), Pp[1] + t * (Q[1] - Pp[1])))
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    if len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


comps = []


def rec(d, poly, idx):
    if d == len(FUN):
        comps.append((poly, tuple(idx)))
        return
    fa, fb = FUN[d]
    vals = [fa * p[0] + fb * p[1] for p in poly]
    mn, mx = min(vals), max(vals)
    for n in range(fl(mn - HI), fl(mx - LO) + 1):
        p2 = clip(poly, fa, fb, n + HI)
        if p2:
            p2 = clip(p2, -fa, -fb, -(n + LO))
        if p2:
            rec(d + 1, p2, idx + [n])


for a in range(P):
    for b in range(P):
        sq = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
        rec(0, sq, [a, -b - 1])


def kappa(idx):
    """max t s.t. f(c) - n >= t and n + 1 - f(c) >= t for all functionals (with the component's indices)"""
    cons = []      # a x + b y + t <= c
    for (fa, fb), n in zip(ALLF, idx):
        cons.append((-fa, -fb, -n))          # t <= f - n
        cons.append((fa, fb, n + 1))         # t <= n + 1 - f
    best = None
    for (c1, c2, c3) in combinations(cons, 3):
        A = [[c[0], c[1], Fr(1)] for c in (c1, c2, c3)]
        bvec = [c[2] for c in (c1, c2, c3)]
        det = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
               + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
        if det == 0:
            continue
        sol = []
        for col in range(3):
            M = [row[:] for row in A]
            for r in range(3):
                M[r][col] = bvec[r]
            dm = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                  + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
            sol.append(dm / det)
        x, y, t = sol
        if all(c[0] * x + c[1] * y + t <= c[2] for c in cons):
            if best is None or t > best:
                best = t
    return best


def is_type(poly, idx):
    # does the component contain a type point (P h or (P/3)(a+bi)) modulo P Z[i]?  test the indices of the point
    pts = [(Fr(P, 2), Fr(P, 2))] + [(Fr(P * a, 3), Fr(P * b, 3)) for a in (1, 2) for b in (1, 2)]
    for T in pts:
        for u in range(-2, 3):
            for v in range(-2, 3):
                c = (T[0] + P * u, T[1] + P * v)
                if tuple(fl(f[0] * c[0] + f[1] * c[1]) for f in ALLF) == idx:
                    return True
    return False


others = [(kappa(idx), idx) for poly, idx in comps if not is_type(poly, idx)]
ntype = len(comps) - len(others)
from collections import Counter
print(f"{name} theta={theta}: {len(comps)} components mod {P}: {ntype} of type c/q, {len(others)} others; "
      f"largest least margins of the others: {dict(Counter(str(k) for k, _ in others))}")
