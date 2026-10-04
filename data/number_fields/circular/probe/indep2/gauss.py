# Exact Gaussian-rational helpers (independent referee code).
from fractions import Fraction as Fr
import math

class G:
    __slots__ = ("re", "im")
    def __init__(self, re, im=0):
        self.re = Fr(re); self.im = Fr(im)
    def __add__(self, o):
        o = o if isinstance(o, G) else G(o)
        return G(self.re + o.re, self.im + o.im)
    __radd__ = __add__
    def __sub__(self, o):
        o = o if isinstance(o, G) else G(o)
        return G(self.re - o.re, self.im - o.im)
    def __rsub__(self, o):
        return G(o) - self
    def __neg__(self):
        return G(-self.re, -self.im)
    def __mul__(self, o):
        o = o if isinstance(o, G) else G(o)
        return G(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__
    def conj(self):
        return G(self.re, -self.im)
    def norm(self):
        return self.re * self.re + self.im * self.im
    def __truediv__(self, o):
        o = o if isinstance(o, G) else G(o)
        n = o.norm()
        p = self * o.conj()
        return G(p.re / n, p.im / n)
    def __pow__(self, e):
        if e < 0:
            return (G(1) / self) ** (-e)
        r = G(1); b = self
        while e:
            if e & 1: r = r * b
            b = b * b; e >>= 1
        return r
    def __eq__(self, o):
        o = o if isinstance(o, G) else G(o)
        return self.re == o.re and self.im == o.im
    def __hash__(self):
        return hash((self.re, self.im))
    def is_gint(self):
        return self.re.denominator == 1 and self.im.denominator == 1
    def __repr__(self):
        return f"({self.re}{'+' if self.im >= 0 else '-'}{abs(self.im)}i)"

I = G(0, 1)
RHO = G(Fr(3, 5), Fr(4, 5))
H = G(Fr(1, 2), Fr(1, 2))

def frac_center(t):
    """t mod 1 in (-1/2, 1/2]."""
    t = Fr(t)
    r = t - math.floor(t)
    if r > Fr(1, 2):
        r -= 1
    return r

def modZ(z):
    """z mod Z[i], coordinates in (-1/2,1/2]."""
    return G(frac_center(z.re), frac_center(z.im))

def mod1(t):
    t = Fr(t)
    return t - math.floor(t)

def dist_int(t):
    """||t||, distance to nearest integer."""
    r = mod1(t)
    return min(r, 1 - r)
