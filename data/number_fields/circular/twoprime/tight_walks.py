"""The finite facts behind the corollary on finite witnesses (circular values below 4 are attained by finite
subgraphs): exact arithmetic in Q(i).

(1) Family 7, window G(2,1): each of the eight dual certificates of value 2/7 in cert_K2M1_seven.txt has three
    positive multipliers; at its 7-point the three margins are exactly 2/7; and the cancellation of their linear parts
    is a relation sum n_i u_i = 0 with positive integers n_i among the tight rotations u_i = +-gamma_i (gamma_i =
    rho^j sigma^l for Re, -i rho^j sigma^l for Im; + for a lower margin, - for an upper one), whose length sum n_i is
    a multiple of 7.
(2) Type q: for each of the four classes alpha + beta i (alpha, beta in {1, 2}) and each sign s, the rotations gamma
    of G_5 with Re(5 gamma (alpha + beta i)) = s (mod 3) contain three with a positive integer relation.
(3) The Moser spindle has chromatic number 4 and maps to K_{7/2}; with 7 vertices its circular chromatic number is
    7/2."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import lcm
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def power(z, e):
    r = (Q(1), Q(0))
    if e < 0:
        n = z[0] * z[0] + z[1] * z[1]
        z, e = (z[0] / n, -z[1] / n), -e
    for _ in range(e):
        r = mul(r, z)
    return r


RHO, SIGMA = (Q(3, 5), Q(4, 5)), (Q(5, 13), Q(12, 13))


def rotation(j, l, part):
    g = mul(power(RHO, j), power(SIGMA, l))
    return g if part == "Re" else mul((Q(0), Q(-1)), g)      # Im(conj(c) g) = Re(conj(c) (-i g))


def value(c, g):                                            # Re(conj(c) g)
    return c[0] * g[0] + c[1] * g[1]


print("(1) family 7 in the window G(2,1): tight relations from the certificates of value 2/7")
lengths = []
for line in open(os.path.join(HERE, "cert_K2M1_seven.txt")).read().splitlines()[2:]:
    head, label, cert = [s.strip() for s in line.split("|")]
    c = tuple(Q(t) for t in label.split()[1:3])
    us, lams = [], []
    for term in cert.split("cert", 1)[1].split(";"):
        j, l, part, sign, lam = term.split()
        g = rotation(int(j), int(l), part)
        v = value(c, g)
        n = v.numerator // v.denominator
        margin = v - n if sign == "+" else n + 1 - v
        assert margin == Q(2, 7)
        eps = 1 if sign == "+" else -1
        us.append((eps * g[0], eps * g[1]))
        lams.append(Q(lam))
        assert (value(c, us[-1]) - Q(2, 7)).denominator == 1      # the tight rotation has value 2/7 mod 1
    assert len(lams) == 3 and all(x > 0 for x in lams) and sum(lams) == 1
    d = lcm(*(x.denominator for x in lams))
    mult = [int(x * d) for x in lams]
    assert (sum(m * u[0] for m, u in zip(mult, us)), sum(m * u[1] for m, u in zip(mult, us))) == (0, 0)
    assert sum(mult) % 7 == 0
    lengths.append(sum(mult))
    print(f"  7-point {c[0]} + {c[1]}i: multiplicities {mult}, length {sum(mult)}")
print(f"  all 8: positive integer relations among tight rotations; lengths {sorted(set(lengths))}")

print("(2) type q: positive relations among the tight rotations of G_5")
G5 = [(Q(1), Q(0)), (Q(-1), Q(0)), (Q(0), Q(1)), (Q(0), Q(-1))]
G5 += [(Q(a * x, 5), Q(b * y, 5)) for x, y in ((3, 4), (4, 3)) for a in (1, -1) for b in (1, -1)]
assert len(set(G5)) == 12 and all(x * x + y * y == 1 for x, y in G5)
for alpha, beta in product((1, 2), repeat=2):
    for s in (1, -1):
        tight = [g for g in G5 if (int(5 * (g[0] * alpha - g[1] * beta)) - s) % 3 == 0]
        assert len(tight) == 6
        found = None
        for a, b, cc in combinations(tight, 3):
            det = a[0] * b[1] - a[1] * b[0]
            if det == 0:
                continue
            x = (-cc[0] * b[1] + cc[1] * b[0]) / det            # x a + y b + cc = 0
            y = (-a[0] * cc[1] + a[1] * cc[0]) / det
            if x > 0 and y > 0:
                d = lcm(x.denominator, y.denominator)
                found = ((a, b, cc), (int(x * d), int(y * d), d))
                break
        assert found is not None
        (a, b, cc), m = found
        assert all(m[0] * a[t] + m[1] * b[t] + m[2] * cc[t] == 0 for t in (0, 1))
        print(f"  class {alpha}+{beta}i, sign {s:+d}: multiplicities {m} on "
              + ", ".join(f"{str(g[0])}{'+' if g[1] >= 0 else '-'}{str(abs(g[1]))}i" for g in (a, b, cc)))

print("(3) the Moser spindle")
EDGES = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4), (0, 5), (4, 5), (4, 6), (5, 6), (3, 6)]
three = [c for c in product(range(3), repeat=7) if all(c[a] != c[b] for a, b in EDGES)]
col = (0, 2, 4, 6, 3, 5, 1)
assert all(2 <= (col[b] - col[a]) % 7 <= 5 for a, b in EDGES)
print(f"  proper 3-colourings: {len(three)}; a homomorphism to K_7/2: {col}; so chi = 4 and chi_c = 7/2")
