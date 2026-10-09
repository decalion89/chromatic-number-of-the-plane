"""k47.py -- exact arithmetic in K = Q(sqrt47)(i) for unit-distance graphs over Q(sqrt47).

An element is a + b r + c i + d i r with r = sqrt47, i^2 = -1, stored as 4 Fractions (a, b, c, d).
Plane: Q(sqrt47)^2 = K via (x, y) -> x + i y, embedded by r -> +sqrt(47). |t| = 1 iff t * conj(t) = 1, with
conj: i -> -i (r fixed). Unit vectors: t = s / conj(s) for s in K^*."""
from fractions import Fraction as Fr
from math import sqrt, gcd
from functools import reduce

import os
R = int(os.environ.get("QD", "47"))


def mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    # (a + b r + c i + d ir)(e + f r + g i + h ir), r^2 = 47, i^2 = -1, (ir)^2 = -47
    return (a * e + R * b * f - c * g - R * d * h,
            a * f + b * e - c * h - d * g,
            a * g + c * e + R * b * h + R * d * f,
            a * h + d * e + b * g + c * f)


def conj(x):
    a, b, c, d = x
    return (a, b, -c, -d)


def rconj(x):
    a, b, c, d = x
    return (a, -b, c, -d)


def inv(x):
    # 1/x = conj(x) * rconj(x conj(x)) / N, where x conj(x) = p + q r in Q(r), N = p^2 - 47 q^2
    xc = mul(x, conj(x))
    p, q = xc[0], xc[1]
    assert xc[2] == 0 and xc[3] == 0
    N = p * p - R * q * q
    num = mul(conj(x), (p, -q, Fr(0), Fr(0)))
    return tuple(Fr(v) / N for v in num)


def norm1(x):
    return mul(x, conj(x)) == (1, 0, 0, 0)


def val(x):
    a, b, c, d = (float(v) for v in x)
    s = sqrt(R)
    return complex(a + b * s, c + d * s)


def unit_from(s):
    s = tuple(Fr(v) for v in s)
    t = mul(s, inv(conj(s)))
    assert norm1(t), (s, t)
    return t
