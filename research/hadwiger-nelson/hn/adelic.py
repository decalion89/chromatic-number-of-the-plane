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


# ---------------------------------------------------------------------------------------------
# Reduction at a non-split place: the field of five_rho7 is 5-colourable.
#
# K = Q(sqrt-3, sqrt-11, sqrt-247), K+ = Q(sqrt33, sqrt741, sqrt2717).  Above 11, K+_v = Q_11(sqrt33)
# (ramified, uniformiser pi = sqrt33, residue field F_11) and K_v = K+_v(sqrt-3) is unramified over it
# (-3 is a non-residue mod 11): the place does not split, every unit vector of K is a unit of O_v
# of relative norm 1, and it reduces to N1 = {x + y w : x^2 + 3 y^2 = 1} in F_121 = F_11(w), w^2 = -3.
# z -> (z mod pi O_v) sends a unit step to a step by an element of N1, so any colouring of the
# finite-plane graph Cay(F_121, N1) (chromatic number 5) colours the unit-distance graph on K.

def _sqrt_padic(n, p, N):
    """s with s^2 = n mod p^N (n a nonzero square mod p, p odd), s mod p the smaller root."""
    r = next(x for x in range(1, p) if (x * x - n) % p == 0)
    s, q = r, p
    for _ in range(N - 1):
        q2 = q * p
        s = (s - (s * s - n) * pow(2 * s, -1, q2)) % q2
        q = q2
    return s


def units_digit(r, p):
    """(r - frac_p(r)) mod p: the p-adic digit of r at p^0."""
    r = Fr(r)
    n, d = r.numerator, r.denominator
    k = 0
    while d % p == 0:
        d //= p
        k += 1
    m = p ** (k + 1)
    return ((n * pow(d, -1, m)) % m) // (p ** k)


_P11 = 11
_S741 = _sqrt_padic(741, 11, 60)


def field247_coordinates(x, y):
    """(a0, a1, a2, a3, b0, b1, b2, b3) with x = a0 + a1 sqrt33 + a2 sqrt741 + a3 sqrt2717 and
    y = b0 sqrt3 + b1 sqrt11 + b2 sqrt247 + b3 sqrt8151, i.e. x + iy in Q(sqrt-3, sqrt-11, sqrt-247).
    A point drawn in i*K is turned by -90 degrees first."""
    F = x.field
    X = dict(zip(F._prod, x.c))
    Y = dict(zip(F._prod, y.c))
    re_, im_ = (1, 33, 741, 2717), (3, 11, 247, 8151)
    if all(X.get(q, 0) == 0 for q in re_) and all(Y.get(q, 0) == 0 for q in im_) and (
            any(X.get(q, 0) for q in im_) or any(Y.get(q, 0) for q in re_)):
        X, Y = Y, {q: -v for q, v in X.items()}
    if any(X.get(q, 0) for q in im_) or any(Y.get(q, 0) for q in re_):
        raise ValueError("not a point of Q(sqrt-3, sqrt-11, sqrt-247) (nor of i times it)")
    return tuple(Fr(X.get(q, 0)) for q in re_) + tuple(Fr(Y.get(q, 0)) for q in im_)


def reduce11(p, place=1):
    """The residue (x0, y0) in F_11^2 of the point p of K at a place above 11 (see above):
    z = alpha + beta sqrt-3 with alpha = a + b pi, beta = c + d pi (pi = sqrt33), and
    (x0, y0) = (units digit of a, units digit of c).  place = +-1 picks sqrt741 in Z_11."""
    a0, a1, a2, a3, b0, b1, b2, b3 = field247_coordinates(p.x, p.y)
    s = _S741 if place == 1 else (11 ** 60 - _S741)
    # sqrt2717 = pi s / 3, sqrt-11 = pi sqrt-3 / 3, sqrt-247 = sqrt-3 s / 3, sqrt-8151 = sqrt-3 pi s / 3
    a = a0 + a2 * s
    c = b0 + b2 * s / 3
    return units_digit(a, 11), units_digit(c, 11)


def finite_plane_11_colouring():
    """A proper 5-colouring of Cay(F_121, N1), N1 = {x + y w : x^2 + 3 y^2 = 1 mod 11}, as a dict."""
    from pysat.solvers import Solver
    q = 11
    N1 = [(a, b) for a in range(q) for b in range(q) if (a * a + 3 * b * b - 1) % q == 0]
    V = [(a, b) for a in range(q) for b in range(q)]
    idx = {v: i for i, v in enumerate(V)}
    k = 5
    var = lambda v, c: v * k + c + 1
    with Solver(name="cadical195") as s:
        for i in range(len(V)):
            s.add_clause([var(i, c) for c in range(k)])
        for (a, b) in V:
            for (x, y) in N1:
                j = idx[((a + x) % q, (b + y) % q)]
                i = idx[(a, b)]
                if i < j:
                    for c in range(k):
                        s.add_clause([-var(i, c), -var(j, c)])
        assert s.solve()
        m = set(l for l in s.get_model() if l > 0)
    return {v: next(c for c in range(k) if var(i, c) in m) for i, v in enumerate(V)}, N1
