"""All unit vectors ((a + b r) / D, (c + e r) / D) of Q(r)^2, r = sqrt(d), with integers a, b, c, e:
a^2 + c^2 + d (b^2 + e^2) = D^2 and a b + c e = 0. For (b, e) != (0, 0), a b + c e = 0 forces
(a, c) = t (e, -b) / g with g = gcd(b, e), and then t^2 = rest g^2 / (b^2 + e^2)."""
from math import gcd, isqrt


def units_fast(d, D):
    out = set()
    D2 = D * D
    for a in range(-D, D + 1):                       # b = e = 0: rational unit vectors
        c2 = D2 - a * a
        c = isqrt(c2)
        if c * c == c2:
            out.add((a, 0, c, 0)); out.add((a, 0, -c, 0))
    B = isqrt(D2 // d)
    for b in range(-B, B + 1):
        for e in range(-B, B + 1):
            if b == 0 and e == 0:
                continue
            s = b * b + e * e
            rest = D2 - d * s
            if rest < 0:
                continue
            g = gcd(b, e)
            num = rest * g * g
            if num % s:
                continue
            t2 = num // s
            t = isqrt(t2)
            if t * t != t2:
                continue
            for tt in ({t, -t} if t else {0}):
                a, c = tt * e // g, -tt * b // g
                out.add((a, b, c, e))
    return sorted(out)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    import os
    for d, D in [(11, 30), (23, 120), (47, 240), (59, 210), (71, 120), (35, 450), (131, 330)]:
        os.environ["QD"] = str(d)
        u = units_fast(d, D)
        # brute force check against the definition
        bf = set()
        D2 = D * D
        B = isqrt(D2 // d)
        for b in range(-B, B + 1):
            for e in range(-B, B + 1):
                rest = D2 - d * (b * b + e * e)
                if rest < 0:
                    continue
                A = isqrt(rest)
                for a in range(-A, A + 1):
                    c2 = rest - a * a
                    c = isqrt(c2)
                    if c * c != c2:
                        continue
                    for cc in ({c, -c} if c else {0}):
                        if a * b + cc * e == 0:
                            bf.add((a, b, cc, e))
        print(d, D, len(u), len(bf), set(u) == bf)
