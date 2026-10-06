"""check_theta.py d D THETA.json: exact check (integers and fractions only).

Lists ALL unit vectors u = ((a + b sqrt d)/D, (c + e sqrt d)/D) of Q(sqrt d)^2 with denominator D
(a^2 + d b^2 + c^2 + d e^2 = D^2 and ab + ce = 0), by brute force over b, e and then a, c, and computes, for the
character  x = ((a + b sqrt d)/D, (c + e sqrt d)/D) -> theta_1 a + theta_2 b + theta_3 c + theta_4 e  (mod 1)
of the group (1/D)Z^4 of such points, the least margin m = min_u ||xi(u)||.  If m > 0, every graph on points of
this group whose edges are unit vectors with denominator D maps to K_{p/q} for every rational p/q with q/p <= m
(Theorem W+, 'if' direction), so its circular chromatic number is at most 1/m."""
import sys, json
from fractions import Fraction as Fr
from math import isqrt


def units(d, D):
    U = set()
    bm = isqrt(D * D // d)
    for b in range(-bm, bm + 1):
        for e in range(-bm, bm + 1):
            r = D * D - d * (b * b + e * e)
            if r < 0:
                continue
            for a in range(-isqrt(r), isqrt(r) + 1):
                c2 = r - a * a
                c = isqrt(c2)
                if c * c != c2:
                    continue
                for cc in {c, -c}:
                    if a * b + cc * e == 0:
                        U.add((a, b, cc, e))
    U = sorted(U)
    assert all(a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0 for a, b, c, e in U)
    return U


def margin(t):
    t = t - (t.numerator // t.denominator)
    return min(t, 1 - t)


def least_margin(d, D, theta):
    th = [Fr(x) for x in theta]
    U = units(d, D)
    return len(U), min(margin(th[0] * a + th[1] * b + th[2] * c + th[3] * e) for a, b, c, e in U)


if __name__ == "__main__":
    d, D = int(sys.argv[1]), int(sys.argv[2])
    theta = json.load(open(sys.argv[3]))["theta"]
    n, m = least_margin(d, D, theta)
    rel = "  > 1/4" if m > Fr(1, 4) else ("  = 1/4" if m == Fr(1, 4) else "")
    print(f"d={d} D={D}: {n} unit vectors (both signs); theta={theta}; least margin {m} = {float(m):.6f}{rel}")
