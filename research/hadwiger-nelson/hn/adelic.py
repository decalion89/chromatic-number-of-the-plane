"""Characters of a number field seen through one place, and the 2-adic 4-colouring of the Moser field.

A circular colouring c(x) = floor(k frac(phi(x))) only needs phi to be additive.  One source of such
phi: complete the field at a place v and take x -> frac_p(L(x)), L a Q_p-linear form on K_v.  If v is
a place of K+ that does not split in the CM field K, every unit vector u of K (u ubar = 1) is a
norm-one element of the local field, and those form a compact group; checking the colouring on that
group checks it on every unit vector of K at once.

For the Moser field K = Q(sqrt-3, sqrt-11): 33 = 1 mod 8, so sqrt33 is a 2-adic integer s and K sits
inside Q_2(omega), omega = (1 + sqrt-3)/2, the unramified quadratic extension of Q_2, where 2 does
not split.  The norm-one units of Z_2[omega] fall into six classes mod 4 (the sixth roots of unity),
and phi(alpha + beta omega) = frac_2((alpha + 2 beta)/4) is 1/4, 1/2 or 3/4 on all of them.
"""
from fractions import Fraction as Fr

PRECISION = 200


def sqrt_2adic(n, N=PRECISION):
    """s with s^2 = n mod 2^N and s = 1 mod 4 (needs n = 1 mod 8)."""
    assert n % 8 == 1
    s = 1
    for k in range(3, N):
        if (s * s - n) % (1 << (k + 1)):
            s += 1 << (k - 1)
    assert (s * s - n) % (1 << N) == 0
    return s if s % 4 == 1 else (1 << N) - s


def frac2(r):
    """The 2-adic fractional part of a rational number, in [0, 1) with a power-of-2 denominator."""
    r = Fr(r)
    n, d = r.numerator, r.denominator
    k = 0
    while d % 2 == 0:
        d //= 2
        k += 1
    if k == 0:
        return Fr(0)
    return Fr((n * pow(d, -1, 1 << k)) % (1 << k), 1 << k)


_S = sqrt_2adic(33)


def moser_coordinates(x, y):
    """(a, b, c, d) with x + iy = a + b sqrt33 + c sqrt-3 + d sqrt-11, for x, y in Q(sqrt3, sqrt11).

    A configuration drawn in i*K (as Exoo and Ismailescu's is) is turned by -90 degrees first; a
    quarter turn keeps every distance.  Raises ValueError for a point outside K and i*K."""
    F = x.field
    X = dict(zip(F._prod, x.c))
    Y = dict(zip(F._prod, y.c))
    if all(X.get(q, 0) == 0 for q in (1, 33)) and all(Y.get(q, 0) == 0 for q in (3, 11)) and any(
            X.get(q, 0) or Y.get(q2, 0) for q, q2 in ((3, 1), (11, 33))):
        X, Y = Y, {q: -v for q, v in X.items()}
    if any(X.get(q, 0) for q in (3, 11)) or any(Y.get(q, 0) for q in (1, 33)):
        raise ValueError("not a point of the Moser field (nor of i times it)")
    return Fr(X.get(1, 0)), Fr(X.get(33, 0)), Fr(Y.get(3, 0)), Fr(Y.get(11, 0))


def moser_phi(p, place=1):
    """frac_2((alpha + 2 beta)/4) for the image alpha + beta omega of p in Q_2(omega).

    place = +1 or -1 picks the square root of 33 (the two places of Q(sqrt33) above 2); both work."""
    a, b, c, d = moser_coordinates(p.x, p.y)
    s = place * _S
    # sqrt-3 = 2 omega - 1 and sqrt-11 = (sqrt33 / 3) sqrt-3
    alpha = a + b * s - c - d * s / 3
    beta = 2 * c + 2 * d * s / 3
    return frac2((alpha + 2 * beta) / 4)


def moser_colour(p, place=1):
    """The colour 0..3 of a point of the Moser field (or of i times it)."""
    return int(4 * moser_phi(p, place))
