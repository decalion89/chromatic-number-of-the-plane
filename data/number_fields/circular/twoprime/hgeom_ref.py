"""Referee's own exact geometry (independent of the probe code).

Convex polygons are kept in H-representation (list of closed half-planes a*x + b*y <= c) and their vertices are found
by intersecting every pair of boundary lines and keeping the feasible intersection points (vertex enumeration), with
redundant half-planes pruned afterwards.  Everything is exact (fractions.Fraction).

Gaussian rationals are pairs (re, im).
"""
from fractions import Fraction as Fr
from math import floor, ceil


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gsub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def gscale(t, a):
    return (t * a[0], t * a[1])


RHO = (Fr(3, 5), Fr(4, 5))


def rho_pow(j):
    z = (Fr(1), Fr(0))
    base = RHO if j >= 0 else gconj(RHO)
    for _ in range(abs(j)):
        z = gmul(z, base)
    return z


def funcs(j):
    """the two linear functionals c -> Re(conj(c) rho^j), Im(conj(c) rho^j), as coefficient pairs"""
    A, B = rho_pow(j)
    # conj(c) rho^j = (c1 - i c2)(A + iB) = (A c1 + B c2) + i (B c1 - A c2)
    return [(A, B), (B, -A)]


def ev(f, p):
    return f[0] * p[0] + f[1] * p[1]


def fl(q):
    return q.numerator // q.denominator


def cl(q):
    return -((-q.numerator) // q.denominator)


def vertices(H):
    """vertices of the bounded convex polygon {p : a p1 + b p2 <= c for (a,b,c) in H} (empty list if empty)"""
    pts = set()
    n = len(H)
    for i in range(n):
        a1, b1, c1 = H[i]
        for j in range(i + 1, n):
            a2, b2, c2 = H[j]
            det = a1 * b2 - a2 * b1
            if det == 0:
                continue
            x = (c1 * b2 - c2 * b1) / det
            y = (a1 * c2 - a2 * c1) / det
            if (x, y) in pts:
                continue
            ok = True
            for (a, b, c) in H:
                if a * x + b * y > c:
                    ok = False
                    break
            if ok:
                pts.add((x, y))
    return sorted(pts)


def prune(H, V):
    """keep only the half-planes that are tight at some vertex (valid for nonempty bounded polygons)"""
    out = []
    seen = set()
    for h in H:
        a, b, c = h
        if any(a * v[0] + b * v[1] == c for v in V):
            # normalise to avoid duplicates
            key = h
            if key not in seen:
                seen.add(key)
                out.append(h)
    return out


class Poly:
    __slots__ = ("H", "V")

    def __init__(self, H, V=None):
        if V is None:
            V = vertices(H)
        self.H = H
        self.V = V
        if V:
            self.H = prune(H, V)

    def empty(self):
        return not self.V

    def translate(self, t):
        H = [(a, b, c + a * t[0] + b * t[1]) for (a, b, c) in self.H]
        V = [(v[0] + t[0], v[1] + t[1]) for v in self.V]
        return Poly(H, V)

    def strip_split(self, f, lo, hi):
        """pieces of self on which f(p) lies in [n + lo, n + hi] for an integer n (lo < hi, hi - lo < 1)"""
        vals = [ev(f, v) for v in self.V]
        mn, mx = min(vals), max(vals)
        out = []
        for n in range(cl(mn - hi), fl(mx - lo) + 1):
            H = self.H + [(f[0], f[1], n + hi), (-f[0], -f[1], -(n + lo))]
            P = Poly(H)
            if not P.empty():
                out.append(P)
        return out

    def ordered(self):
        V = self.V
        if len(V) <= 2:
            return list(V)
        cx = sum(v[0] for v in V) / len(V)
        cy = sum(v[1] for v in V) / len(V)
        # exact angular sort via half-plane + cross product
        def half(v):
            dx, dy = v[0] - cx, v[1] - cy
            return 0 if (dy > 0 or (dy == 0 and dx > 0)) else 1
        from functools import cmp_to_key
        def cmp(u, v):
            hu, hv = half(u), half(v)
            if hu != hv:
                return hu - hv
            cr = (u[0] - cx) * (v[1] - cy) - (u[1] - cy) * (v[0] - cx)
            return -1 if cr > 0 else (1 if cr < 0 else 0)
        return sorted(V, key=cmp_to_key(cmp))

    def area(self):
        V = self.ordered()
        if len(V) < 3:
            return Fr(0)
        s = Fr(0)
        for i in range(len(V)):
            x1, y1 = V[i]
            x2, y2 = V[(i + 1) % len(V)]
            s += x1 * y2 - x2 * y1
        return abs(s) / 2

    def centroid_v(self):
        V = self.V
        return (sum(v[0] for v in V) / len(V), sum(v[1] for v in V) / len(V))
