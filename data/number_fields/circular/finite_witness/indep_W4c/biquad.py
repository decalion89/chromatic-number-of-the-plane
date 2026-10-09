"""Exact arithmetic in the biquadratic field Q(sqrt(A), sqrt(B)).

An element is a 4-tuple of Fractions (x0, x1, x2, x3) meaning
x0 + x1*sa + x2*sb + x3*sab  with sa^2 = A, sb^2 = B, sab = sa*sb.

The products are
  sa*sa = A, sb*sb = B, sab*sab = A*B,
  sa*sb = sab, sa*sab = A*sb, sb*sab = B*sa.

Integer-only helpers (int_mul, int_sq_norm) work on integer 4-tuples and are
used for the bulk exact checks: everything is plain Python int arithmetic
(arbitrary precision), so no rounding can occur anywhere.
"""
from fractions import Fraction


class Field:
    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.AB = A * B

    # ---- generic (Fraction or int) arithmetic on 4-tuples -------------------
    def mul(self, x, y):
        A, B, AB = self.A, self.B, self.AB
        x0, x1, x2, x3 = x
        y0, y1, y2, y3 = y
        return (
            x0 * y0 + A * x1 * y1 + B * x2 * y2 + AB * x3 * y3,
            x0 * y1 + x1 * y0 + B * (x2 * y3 + x3 * y2),
            x0 * y2 + x2 * y0 + A * (x1 * y3 + x3 * y1),
            x0 * y3 + x3 * y0 + x1 * y2 + x2 * y1,
        )

    @staticmethod
    def add(x, y):
        return tuple(a + b for a, b in zip(x, y))

    @staticmethod
    def sub(x, y):
        return tuple(a - b for a, b in zip(x, y))

    @staticmethod
    def neg(x):
        return tuple(-a for a in x)

    def const(self, q):
        return (Fraction(q), Fraction(0), Fraction(0), Fraction(0))

    def to_float(self, x):
        import math
        sa, sb = math.sqrt(self.A), math.sqrt(self.B)
        return float(x[0]) + float(x[1]) * sa + float(x[2]) * sb + float(x[3]) * sa * sb

    # ---- vectors in F^2 viewed as complex numbers (re, im) -------------------
    def cmul(self, u, v):
        """(u0 + i u1)(v0 + i v1) for u, v in F^2."""
        a, b = u
        c, d = v
        return (self.sub(self.mul(a, c), self.mul(b, d)),
                self.add(self.mul(a, d), self.mul(b, c)))

    def sq_norm(self, u):
        a, b = u
        return self.add(self.mul(a, a), self.mul(b, b))


def int_sq_norm_8(F, d):
    """Exact squared Euclidean norm of the vector with integer numerators
    d = (a0..a3, b0..b3) (common denominator left out).  Returns the integer
    4-tuple of  (a0+a1 sa+a2 sb+a3 sab)^2 + (b0+...)^2."""
    A, B, AB = F.A, F.B, F.AB
    a0, a1, a2, a3, b0, b1, b2, b3 = d
    return (
        a0 * a0 + A * a1 * a1 + B * a2 * a2 + AB * a3 * a3
        + b0 * b0 + A * b1 * b1 + B * b2 * b2 + AB * b3 * b3,
        2 * (a0 * a1 + B * a2 * a3 + b0 * b1 + B * b2 * b3),
        2 * (a0 * a2 + A * a1 * a3 + b0 * b2 + A * b1 * b3),
        2 * (a0 * a3 + a1 * a2 + b0 * b3 + b1 * b2),
    )


def vec_to_int8(u, D):
    """Convert a vector of two Fraction 4-tuples to an integer 8-tuple of
    numerators over D; raises if not exactly representable."""
    out = []
    for comp in u:
        for c in comp:
            q = Fraction(c) * D
            if q.denominator != 1:
                raise ValueError("not representable over denominator %d: %s" % (D, c))
            out.append(int(q))
    return tuple(out)


def int8_to_vec(t, D):
    return (tuple(Fraction(t[i], D) for i in range(4)),
            tuple(Fraction(t[i], D) for i in range(4, 8)))
