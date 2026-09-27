"""Square roots in a multiquadratic field, done by characters not by search.

An element t of Q(sqrt d1, ..., sqrt dk) is a square exactly when some y in
the field has y^2 = t.  Numerically, the conjugates of y must be plus or
minus the square roots of the conjugates of t -- and the sign pattern is not
free: it has to be a character of the Galois group (Z/2)^k, so there are 2^k
patterns to try, not 2^(2^k - 1).  Getting that wrong turns sixteen linear
solves into thirty-two thousand.

Each candidate pattern gives the conjugates of y, one linear solve gives its
rational coordinates, and the answer is accepted only when squaring it
returns t exactly -- so the numerics only propose, and the field decides.
"""
import itertools
import numpy as np
from fractions import Fraction as Fr


class Roots:
    def __init__(self, field):
        self.K = field
        self.gens = field.gens
        self.dim = field.dim
        k = len(self.gens)
        self.embs = list(itertools.product(*[(1, -1)] * k))
        B = np.zeros((len(self.embs), self.dim))
        for si, sg in enumerate(self.embs):
            for m in range(self.dim):
                v = 1.0
                for i, g in enumerate(self.gens):
                    if m >> i & 1:
                        v *= sg[i] * (g ** .5)
                B[si, m] = v
        self.B = B
        self.Binv = np.linalg.inv(B)
        # the 2^k characters: chi_S(sigma) = product of sigma's signs over S
        self.chars = []
        for S in range(self.dim):
            col = np.ones(len(self.embs))
            for si, sg in enumerate(self.embs):
                v = 1
                for i in range(k):
                    if S >> i & 1:
                        v *= sg[i]
                col[si] = v
            self.chars.append(col)

    def value(self, e):
        return self.B @ np.array([float(x) for x in e.c])

    def sqrt(self, e, den=1 << 16, tol=1e-7):
        vals = self.value(e)
        if (vals < -1e-12).any():
            return None
        root = np.sqrt(np.maximum(vals, 0.0))
        for chi in self.chars:
            coef = self.Binv @ (root * chi)
            cand = []
            ok = True
            for x in coef:
                f = Fr(x).limit_denominator(den)
                if abs(float(f) - x) > tol:
                    ok = False
                    break
                cand.append(f)
            if not ok:
                continue
            y = self.K.zero()
            for m, f in enumerate(cand):
                if not f:
                    continue
                t = self.K.rational(f)
                for i, g in enumerate(self.gens):
                    if m >> i & 1:
                        t = t * self.K.sqrt(g)
                y = y + t
            if y * y == e:
                return y
        return None
