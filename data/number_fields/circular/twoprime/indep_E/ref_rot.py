"""Referee E: shared tiny helpers (rotations only).  Exact integer arithmetic.

A rational rotation gamma with 65*gamma in Z[i] is stored as the integer pair (p, q),
gamma = (p + q i)/65.  For c = x + i y:  conj(c)*gamma = ((p x + q y) + i (q x - p y))/65,
so Re(conj(c) gamma) = (p x + q y)/65 and Im(conj(c) gamma) = Re(conj(c) * (-i gamma)).
"""


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def rotations_from_generators(M5=1, M13=1):
    """All i^a rho^j sigma^l, |j|<=M5, |l|<=M13, as integer numerators over D=5^M5*13^M13.
    rho=(2+i)^2/5, sigma=(3+2i)^2/13: D*rho^j*sigma^l = (2+i)^(M5+j)(2-i)^(M5-j)(3+2i)^(M13+l)(3-2i)^(M13-l)."""
    D = 5 ** M5 * 13 ** M13
    out = {}
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    for j in range(-M5, M5 + 1):
        for l in range(-M13, M13 + 1):
            z = (1, 0)
            for _ in range(M5 + j):
                z = gmul(z, (2, 1))
            for _ in range(M5 - j):
                z = gmul(z, (2, -1))
            for _ in range(M13 + l):
                z = gmul(z, (3, 2))
            for _ in range(M13 - l):
                z = gmul(z, (3, -2))
            assert z[0] ** 2 + z[1] ** 2 == D * D
            for a, u in enumerate(units):
                out[(a, j, l)] = gmul(u, z)
    return D, out


def all_norm_solutions(D):
    return sorted((p, q) for p in range(-D, D + 1) for q in range(-D, D + 1) if p * p + q * q == D * D)


def eighteen_functionals():
    """One representative (p,q) per +-pair among the 36 rotations with denominator dividing 65,
    built from generators; checked to be all solutions of p^2+q^2=65^2."""
    D, rots = rotations_from_generators(1, 1)
    vals = sorted(set(rots.values()))
    assert len(vals) == 36, len(vals)
    assert vals == all_norm_solutions(65)
    reps = sorted(set(v if (v[0] > 0 or (v[0] == 0 and v[1] > 0)) else (-v[0], -v[1]) for v in vals))
    assert len(reps) == 18
    return reps


if __name__ == "__main__":
    print(eighteen_functionals())
