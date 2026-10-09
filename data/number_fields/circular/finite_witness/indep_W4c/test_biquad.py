#!/usr/bin/env python3
"""Self-test of the exact field arithmetic against floating point, and of
basic identities (sa^2 = A, sa*sb = sab, ...)."""
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from biquad import Field, int_sq_norm_8  # noqa: E402

random.seed(12345)
for (A, B) in ((3, 11), (2, 3)):
    F = Field(A, B)
    sa = (0, 1, 0, 0)
    sb = (0, 0, 1, 0)
    sab = (0, 0, 0, 1)
    assert F.mul(sa, sa) == (A, 0, 0, 0)
    assert F.mul(sb, sb) == (B, 0, 0, 0)
    assert F.mul(sab, sab) == (A * B, 0, 0, 0)
    assert F.mul(sa, sb) == (0, 0, 0, 1)
    assert F.mul(sa, sab) == (0, 0, A, 0)
    assert F.mul(sb, sab) == (0, B, 0, 0)
    worst = 0.0
    for _ in range(20000):
        x = tuple(Fraction(random.randint(-50, 50), random.randint(1, 9)) for _ in range(4))
        y = tuple(Fraction(random.randint(-50, 50), random.randint(1, 9)) for _ in range(4))
        z = F.mul(x, y)
        assert F.mul(y, x) == z
        fz = F.to_float(z)
        fx, fy = F.to_float(x), F.to_float(y)
        worst = max(worst, abs(fz - fx * fy) / (1 + abs(fx * fy)))
        d = tuple(random.randint(-300, 300) for _ in range(8))
        n8 = int_sq_norm_8(F, d)
        a = d[:4]
        b = d[4:]
        assert n8 == F.add(F.mul(a, a), F.mul(b, b))
    print("A=%d B=%d: identities ok, 20000 random products agree with floats (max rel err %.2e)" % (A, B, worst))
print("biquad self-test PASSED")
