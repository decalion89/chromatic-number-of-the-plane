"""Exact kappa(S) = max_xi min_{s in S} ||xi(s)|| for finite S in Gamma = Z^r x Z/m (r = 1 or 2).

Characters: xi(x, y) = alpha . x + b*y/m (mod 1), alpha in R^r / Z^r, b in Z/m.
For fixed b the maximum of g(alpha) = min_s ||alpha . x_s + c_s|| is attained at a vertex of an LP
  max t  s.t.  alpha . x_s + c_s - n_s >= t,  n_s + 1 - alpha . x_s - c_s >= t,
so it suffices to enumerate solutions of r+1 equations  alpha . x_s - sigma t = n - c_s
(sigma = +-1, n integer) and evaluate g exactly there (Fractions).  The free parts of S span R^r
(S generates Gamma), so the LP is pointed and the optimum is at such a vertex.

Also: the tight set T(xi) = {s in S : xi(s) = kappa mod 1} and an exact test for a positive integer
relation sum n_t t = 0 (n_t >= 0, not all 0) in Gamma (torsion included).
"""
from fractions import Fraction as Fr
from itertools import combinations, product
import math


def frac_part(v):
    return v - math.floor(v)


def dist_int(v):
    f = frac_part(v)
    return min(f, 1 - f)


class Group:
    def __init__(self, r, m):
        self.r, self.m = r, m  # Gamma = Z^r x Z/m  (m = 1: no torsion)

    def norm(self, s):
        return tuple(s[:self.r]) + ((s[self.r] % self.m),) if self.m > 1 else tuple(s[:self.r]) + (0,)

    def neg(self, s):
        return self.norm(tuple(-v for v in s[:self.r]) + (-s[self.r],))

    def add(self, s, t):
        return self.norm(tuple(a + b for a, b in zip(s[:self.r], t[:self.r])) + (s[self.r] + t[self.r],))

    def xi(self, alpha, b, s):
        val = sum(Fr(a) * x for a, x in zip(alpha, s[:self.r])) + Fr(b * s[self.r], self.m)
        return frac_part(val)


def symmetrize(G, S_half):
    S = set()
    for s in S_half:
        s = G.norm(s)
        S.add(s)
        S.add(G.neg(s))
    return sorted(S)


def half(G, S):
    out, seen = [], set()
    for s in S:
        if s in seen:
            continue
        seen.add(s)
        seen.add(G.neg(s))
        out.append(s)
    return out


def solve(Mat, rhs):
    """Exact Gaussian elimination; returns None if singular."""
    n = len(Mat)
    A = [list(map(Fr, row)) + [Fr(v)] for row, v in zip(Mat, rhs)]
    for col in range(n):
        piv = next((i for i in range(col, n) if A[i][col] != 0), None)
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        for i in range(n):
            if i != col and A[i][col] != 0:
                f = A[i][col] / A[col][col]
                A[i] = [a - f * c for a, c in zip(A[i], A[col])]
    return [A[i][n] / A[i][i] for i in range(n)]


def gval(G, Sh, alpha, b):
    return min(dist_int(G.xi(alpha, b, s)) for s in Sh)


def kappa(G, S):
    """Return (kappa, list of optimal vertices (alpha mod 1, b))."""
    Sh = half(G, S)
    r = G.r
    best = Fr(-1)
    opt = set()
    for b in range(G.m):
        cs = {s: Fr(b * s[r], G.m) for s in Sh}
        cons = []
        for s in Sh:
            bound = sum(abs(x) for x in s[:r]) + 2
            for sigma in (1, -1):
                for n in range(-bound, bound + 1):
                    cons.append((s, sigma, n))
        for combo in combinations(cons, r + 1):
            Mat = [list(s[:r]) + [-sigma] for (s, sigma, n) in combo]
            rhs = [n - cs[s] for (s, sigma, n) in combo]
            sol = solve(Mat, rhs)
            if sol is None:
                continue
            alpha = tuple(frac_part(a) for a in sol[:r])
            t = sol[r]
            if t < 0 or t > Fr(1, 2):
                continue
            g = gval(G, Sh, alpha, b)
            if g > best:
                best = g
                opt = {(alpha, b)}
            elif g == best:
                opt.add((alpha, b))
    return best, sorted(opt)


def tight_set(G, S, alpha, b, k):
    return [s for s in S if G.xi(alpha, b, s) == k]


def positive_relation(G, T):
    """Exact test: is there n >= 0 (not all 0) with sum n_t t = 0 in Gamma = Z^r x Z/m?
    Equivalent to 0 in conv of the free parts (then multiply an integral relation by m).
    Returns a relation (dict t -> n) or None."""
    r = G.r
    if not T:
        return None
    # torsion element alone
    for t in T:
        if all(v == 0 for v in t[:r]):
            return {t: G.m}
    pts = [tuple(Fr(v) for v in t[:r]) for t in T]
    if r == 1:
        pos = [t for t in T if t[0] > 0]
        neg = [t for t in T if t[0] < 0]
        if pos and neg:
            a, c = pos[0], neg[0]
            l = math.lcm(a[0], -c[0])
            rel = {a: l // a[0], c: l // (-c[0])}
            return scale_to_group(G, rel)
        return None
    if r == 2:
        def cross(p, q):
            return p[0] * q[1] - p[1] * q[0]

        def dot(p, q):
            return p[0] * q[0] + p[1] * q[1]
        n = len(T)
        for i in range(n):
            for j in range(i + 1, n):
                p, q = pts[i], pts[j]
                if cross(p, q) == 0 and dot(p, q) < 0:
                    # n1 p + n2 q = 0
                    # find integers: p = -lam q, lam > 0
                    lam = (-p[0] / q[0]) if q[0] != 0 else (-p[1] / q[1])
                    rel = {T[i]: lam.denominator, T[j]: lam.numerator}
                    return scale_to_group(G, rel)
        for i, j, k in combinations(range(n), 3):
            p, q, w = pts[i], pts[j], pts[k]
            # barycentric: 0 = a p + b q + c w, a,b,c > 0
            d1, d2, d3 = cross(q, w), cross(w, p), cross(p, q)
            # 0 inside triangle iff d1,d2,d3 all same strict sign (then a=d1,b=d2,c=d3 works)
            if (d1 > 0 and d2 > 0 and d3 > 0) or (d1 < 0 and d2 < 0 and d3 < 0):
                coeffs = [abs(d1), abs(d2), abs(d3)]
                l = 1
                for cf in coeffs:
                    l = math.lcm(l, cf.denominator)
                ints = [int(cf * l) for cf in coeffs]
                g = math.gcd(*ints)
                rel = {T[i]: ints[0] // g, T[j]: ints[1] // g, T[k]: ints[2] // g}
                return scale_to_group(G, rel)
        return None
    raise NotImplementedError


def scale_to_group(G, rel):
    """Given a relation that vanishes in the free part, multiply it so that it vanishes in Gamma."""
    r = G.r
    tot = (0,) * r + (0,)
    for t, n in rel.items():
        for _ in range(n):
            tot = G.add(tot, t)
    assert all(v == 0 for v in tot[:r]), (rel, tot)
    if tot[r] % G.m == 0:
        return rel
    # multiply by the order of the torsion part
    o = G.m // math.gcd(G.m, tot[r])
    return {t: n * o for t, n in rel.items()}


def check_relation(G, rel):
    tot = (0,) * G.r + (0,)
    for t, n in rel.items():
        assert n > 0
        for _ in range(n):
            tot = G.add(tot, t)
    return all(v == 0 for v in tot[:G.r]) and tot[G.r] % G.m == 0
