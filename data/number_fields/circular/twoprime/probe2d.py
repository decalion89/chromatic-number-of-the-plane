"""The two-prime probe (exact).

For N = 5^k 13^m let G(k, m) = { u rho^j sigma^l : u in mu_4, |j| <= k, |l| <= m }, rho = (3+4i)/5,
sigma = (3+2i)/(3-2i) = (5+12i)/13.  These are rational rotations (unit vectors of Q(i)) with denominators dividing N,
and the Z-span of G(k, m) is (1/N) Z[i] (the gcd argument of Lemma 4(i) of the paper at each prime).  So a character
of (1/N) Z[i] is a point c in C / N Z[i] (a -> Re(conj(c) a)), and

    S^r(k, m) = { c in C : Re(conj(c) g) in [r, 1 - r] + Z for all g in G(k, m) }      (periodic mod N Z[i])

is a finite union of convex polygons mod N (each component is the intersection of one strip from each family).
Every character of Q(i) whose values on all rational rotations lie in [r, 1-r] restricts to a point of S^r(k, m).

Lifting: S^r(k+1, m) is obtained from S^r(k, m) by translating each component by N(a + bi), 0 <= a, b <= 4, and
cutting by the conditions at j = +-(k+1), |l| <= m; S^r(k, m+1) likewise with 0 <= a, b <= 12 and l = +-(m+1).

Main components: the component of N h (h = (1+i)/2) and those of the four type q points (N/3)(alpha + beta i),
alpha, beta in {1, 2}; membership is tested exactly through the strip indices (as in levels3.py).
"""
from fractions import Fraction as Fr
import sys, time

ONE = Fr(1)
RHO = (Fr(3, 5), Fr(4, 5))
SIG = (Fr(5, 13), Fr(12, 13))


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def conj(a):
    return (a[0], -a[1])


def cpow(z, n):
    w = (ONE, Fr(0))
    base = z if n >= 0 else conj(z)          # |z| = 1, so z^-1 = conj(z)
    for _ in range(abs(n)):
        w = cmul(w, base)
    return w


def gamma(j, l):
    return cmul(cpow(RHO, j), cpow(SIG, l))


def floor(q):
    return q.numerator // q.denominator


def clip(poly, a, b, c, sense):
    """keep sense*(a x + b y - c) >= 0"""
    out = []
    n = len(poly)
    vals = [sense * (a * p[0] + b * p[1] - c) for p in poly]
    if n == 1:
        return list(poly) if vals[0] >= 0 else []
    for i in range(n):
        P, Q = poly[i], poly[(i + 1) % n]
        fp, fq = vals[i], vals[(i + 1) % n]
        if fp >= 0:
            out.append(P)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    while len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


def strip_split(poly, a, b, lo, hi):
    vals = [a * p[0] + b * p[1] for p in poly]
    vmin, vmax = min(vals), max(vals)
    out = []
    for n in range(floor(vmin - hi), floor(vmax - lo) + 1):
        if n + hi < vmin or n + lo > vmax:
            continue
        q = clip(poly, a, b, n + lo, +1)
        if q:
            q = clip(q, a, b, n + hi, -1)
        if q:
            out.append(q)
    return out


def families(gs):
    """for g = A + iB: Re(conj(c) g) = A x + B y, Im(conj(c) g) = B x - A y"""
    out = []
    for (A, B) in gs:
        out.append((A, B))
        out.append((B, -A))
    return out


def centroid(poly):
    n = len(poly)
    return (sum(p[0] for p in poly) / n, sum(p[1] for p in poly) / n)


class Probe:
    def __init__(self, r):
        self.r = Fr(r)
        lo, hi = self.r, 1 - self.r
        self.k = 0
        self.m = 0
        self.comps = [[(lo, lo), (hi, lo), (hi, hi), (lo, hi)]]

    @property
    def N(self):
        return 5 ** self.k * 13 ** self.m

    def lift5(self):
        lo, hi = self.r, 1 - self.r
        N = self.N
        k1 = self.k + 1
        fams = families([gamma(j, l) for j in (k1, -k1) for l in range(-self.m, self.m + 1)])
        new = []
        for C in self.comps:
            for a in range(5):
                for b in range(5):
                    pieces = [[(p[0] + N * a, p[1] + N * b) for p in C]]
                    for (fa, fb) in fams:
                        nxt = []
                        for P in pieces:
                            nxt.extend(strip_split(P, fa, fb, lo, hi))
                        pieces = nxt
                        if not pieces:
                            break
                    new.extend(pieces)
        self.k = k1
        self.comps = new

    def lift13(self):
        lo, hi = self.r, 1 - self.r
        N = self.N
        m1 = self.m + 1
        fams = families([gamma(j, l) for l in (m1, -m1) for j in range(-self.k, self.k + 1)])
        new = []
        for C in self.comps:
            for a in range(13):
                for b in range(13):
                    pieces = [[(p[0] + N * a, p[1] + N * b) for p in C]]
                    for (fa, fb) in fams:
                        nxt = []
                        for P in pieces:
                            nxt.extend(strip_split(P, fa, fb, lo, hi))
                        pieces = nxt
                        if not pieces:
                            break
                    new.extend(pieces)
        self.m = m1
        self.comps = new

    def all_families(self, k=None, m=None):
        k = self.k if k is None else k
        m = self.m if m is None else m
        return families([gamma(j, l) for j in range(-k, k + 1) for l in range(-m, m + 1)])

    def type_points(self, k=None, m=None):
        k = self.k if k is None else k
        m = self.m if m is None else m
        N = 5 ** k * 13 ** m
        return [('C', (Fr(N, 2), Fr(N, 2)))] + [('Q', (Fr(N * a, 3), Fr(N * b, 3))) for a in (1, 2) for b in (1, 2)]

    def label(self, pt, k=None, m=None):
        """'C', 'Q' or 'X': is the point in the same component (at window (k, m)) as a type point?"""
        lo = self.r
        fams = self.all_families(k, m)
        for name, TP in self.type_points(k, m):
            if all(floor(fa * TP[0] + fb * TP[1] - lo) == floor(fa * pt[0] + fb * pt[1] - lo) for (fa, fb) in fams):
                return name
        return 'X'


def kappa_point(pt, k, m):
    best = Fr(1, 2)
    for (fa, fb) in families([gamma(j, l) for j in range(-k, k + 1) for l in range(-m, m + 1)]):
        v = fa * pt[0] + fb * pt[1]
        f = v - floor(v)
        best = min(best, f, 1 - f)
    return best
