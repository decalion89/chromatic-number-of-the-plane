"""Exact arithmetic in F = Q(t_1, ..., t_r) given as a tensor product of simple extensions Q[x]/(f_j) (assumed
linearly disjoint), and in L = F(i) = F (x) Q(i).  Elements of F: dict {exponent tuple: Fraction}; of L: dict with
Gaussian-rational values (pairs of Fractions)."""
from fractions import Fraction as Fr
import itertools


class Field:
    def __init__(self, polys):
        # polys: list of coefficient lists, low degree first, monic, e.g. x^2 - 3 -> [-3, 0, 1]
        self.polys = [[Fr(c) for c in p] for p in polys]
        self.degs = [len(p) - 1 for p in polys]
        self.basis = list(itertools.product(*[range(d) for d in self.degs]))
        self.index = {b: k for k, b in enumerate(self.basis)}
        self.n = len(self.basis)
        # reduction of t_j^m for m < 2 deg_j - 1: express as poly of degree < deg_j
        self.red = []
        for p, d in zip(self.polys, self.degs):
            table = {}
            for m in range(2 * d - 1):
                v = [Fr(0)] * (m + 1); v[m] = Fr(1)
                for k in range(m, d - 1, -1):          # reduce x^k using x^d = -sum p_i x^i
                    c = v[k]
                    if c:
                        v[k] = Fr(0)
                        for i2 in range(d):
                            v[k - d + i2] -= c * p[i2]
                table[m] = v[:d] + [Fr(0)] * (d - len(v[:d]))
            self.red.append(table)

    def vec(self, a):
        return [a.get(b, Fr(0)) for b in self.basis]

    def mul(self, a, b):
        out = {}
        for ea, ca in a.items():
            if not ca: continue
            for eb, cb in b.items():
                if not cb: continue
                # product of monomials: componentwise, each reduced
                parts = []
                for j, (x, y) in enumerate(zip(ea, eb)):
                    parts.append(self.red[j][x + y])
                for combo in itertools.product(*[list(enumerate(pt)) for pt in parts]):
                    coef = ca * cb
                    e = []
                    for (k, c) in combo:
                        coef *= c; e.append(k)
                        if not coef: break
                    if coef:
                        e = tuple(e); out[e] = out.get(e, Fr(0)) + coef
        return {k: v for k, v in out.items() if v}

    def add(self, a, b, s=1):
        out = dict(a)
        for k, v in b.items():
            out[k] = out.get(k, Fr(0)) + s * v
        return {k: v for k, v in out.items() if v}

    def scal(self, c, a):
        return {k: c * v for k, v in a.items() if c * v}

    def one(self):
        return {tuple([0] * len(self.degs)): Fr(1)}

    def const(self, c):
        return {tuple([0] * len(self.degs)): Fr(c)} if c else {}

    def gen(self, j):
        e = [0] * len(self.degs); e[j] = 1
        return {tuple(e): Fr(1)}

    def inv(self, a):
        """solve a * x = 1 by linear algebra over Q"""
        n = self.n
        cols = []
        for b in self.basis:
            cols.append(self.vec(self.mul(a, {b: Fr(1)})))
        # matrix A with A[:, k] = cols[k]; solve A x = e_0
        A = [[cols[k][r] for k in range(n)] + [Fr(1) if r == 0 else Fr(0)] for r in range(n)]
        for c in range(n):
            piv = next(r for r in range(c, n) if A[r][c] != 0)
            A[c], A[piv] = A[piv], A[c]
            pv = A[c][c]
            A[c] = [x / pv for x in A[c]]
            for r in range(n):
                if r != c and A[r][c] != 0:
                    f = A[r][c]
                    A[r] = [x - f * y for x, y in zip(A[r], A[c])]
        return {self.basis[k]: A[k][n] for k in range(n) if A[k][n]}


def unit_from(Fd, p):
    """v = (p + i)^2/(p^2 + 1) in L = F(i), returned as (Re part, Im part) in F: ((p^2 - 1)/(p^2+1), 2p/(p^2+1))"""
    p2 = Fd.mul(p, p)
    den = Fd.inv(Fd.add(p2, Fd.one()))
    X = Fd.mul(Fd.add(p2, Fd.one(), -1), den)
    Y = Fd.mul(Fd.scal(Fr(2), p), den)
    # check X^2 + Y^2 = 1
    assert Fd.add(Fd.add(Fd.mul(X, X), Fd.mul(Y, Y)), Fd.one(), -1) == {}, "not a unit vector"
    return (X, Y)


def lmul(Fd, u, v):
    """product in L of u = (X1, Y1), v = (X2, Y2) meaning X + iY"""
    return (Fd.add(Fd.mul(u[0], v[0]), Fd.mul(u[1], v[1]), -1), Fd.add(Fd.mul(u[0], v[1]), Fd.mul(u[1], v[0])))


def lconj(Fd, u):
    return (u[0], Fd.scal(Fr(-1), u[1]))
