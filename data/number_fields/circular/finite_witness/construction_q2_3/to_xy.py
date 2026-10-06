"""to_xy.py: coordinates of zeta_24^k = cos(15k deg) + i sin(15k deg), k = 0..7, over the basis (1, sqrt2, sqrt3, sqrt6)
of Q(sqrt2, sqrt3), times 4; and the map from power-basis coordinates (over D) to plane coordinates (over 4D)."""
from fractions import Fraction as Fr
# cos and sin of 15k degrees, k = 0..7, times 4, over (1, s2, s3, s6)
COS4 = [(4, 0, 0, 0), (0, 1, 0, 1), (0, 0, 2, 0), (0, 2, 0, 0), (2, 0, 0, 0), (0, -1, 0, 1), (0, 0, 0, 0), (0, 1, 0, -1)]
SIN4 = [(0, 0, 0, 0), (0, -1, 0, 1), (2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 1, 0, 1), (4, 0, 0, 0), (0, 1, 0, 1)]


def xy(c):
    """power-basis integer coordinates c_0..c_7 (over D) -> 8 integers (x over (1,s2,s3,s6), y ...) over 4D"""
    x = [sum(c[k] * COS4[k][t] for k in range(8)) for t in range(4)]
    y = [sum(c[k] * SIN4[k][t] for k in range(8)) for t in range(4)]
    return x + y


def sq(p, q, x0, x1, x2, x3):
    return (x0 * x0 + p * x1 * x1 + q * x2 * x2 + p * q * x3 * x3, 2 * x0 * x1 + 2 * q * x2 * x3,
            2 * x0 * x2 + 2 * p * x1 * x3, 2 * x0 * x3 + 2 * x1 * x2)


if __name__ == "__main__":
    import math, json
    from z24 import is_unit
    # numeric check of the tables
    s2, s3, s6 = math.sqrt(2), math.sqrt(3), math.sqrt(6)
    val = lambda v: (v[0] + v[1] * s2 + v[2] * s3 + v[3] * s6) / 4
    for k in range(8):
        assert abs(val(COS4[k]) - math.cos(math.radians(15 * k))) < 1e-12, k
        assert abs(val(SIN4[k]) - math.sin(math.radians(15 * k))) < 1e-12, k
    C = json.load(open("cfg_s2_2.json"))
    D = C["D"]
    for u in C["units"]:
        v = xy(u)
        a, b = sq(2, 3, *v[:4]), sq(2, 3, *v[4:])
        assert tuple(x + y for x, y in zip(a, b)) == (16 * D * D, 0, 0, 0), u
        assert is_unit(u, D)
    print("tables and units OK:", len(C["units"]), "vectors, plane denominator", 4 * D)
