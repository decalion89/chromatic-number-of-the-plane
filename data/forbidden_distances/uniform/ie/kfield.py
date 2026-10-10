"""Exact arithmetic in K = Q(sqrt3, sqrt11): elements a + b sqrt3 + c sqrt11 + d sqrt33 with rational a, b, c, d."""
from fractions import Fraction as F
class K:
    __slots__ = ('v',)
    def __init__(self, v):
        self.v = tuple(F(x) for x in v)
    def __add__(s, o): return K(tuple(a + b for a, b in zip(s.v, o.v)))
    def __sub__(s, o): return K(tuple(a - b for a, b in zip(s.v, o.v)))
    def __mul__(s, o):
        a, b, c, d = s.v; e, f, g, h = o.v
        # 1, r3, r11, r33 with r3^2 = 3, r11^2 = 11, r33^2 = 33, r3 r11 = r33, r3 r33 = 3 r11, r11 r33 = 11 r3
        return K((a*e + 3*b*f + 11*c*g + 33*d*h,
                  a*f + b*e + 11*c*h + 11*d*g,
                  a*g + c*e + 3*b*h + 3*d*f,
                  a*h + d*e + b*g + c*f))
    def __eq__(s, o): return s.v == o.v
    def __hash__(s): return hash(s.v)
    def __repr__(s): return 'K(%s)' % ', '.join(str(x) for x in s.v)
ONE = K((1, 0, 0, 0))
def dist2(p, q):
    """|p - q|^2 for points given as (re, im) with re, im in K (both real)"""
    dx = p[0] - q[0]; dy = p[1] - q[1]
    return dx * dx + dy * dy
