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


# ---------------------------------------------------------------------------------------------
# The plane Q(sqrt3, sqrt11)^2 is 4-colourable, so its chromatic number is 4 (the Moser spindle).
#
# A point (x, y) with x, y in L = Q(sqrt3, sqrt11) is z = x + iy in K = L(i) = Q(i, sqrt3, sqrt11).
# Above 2, sqrt33 is a 2-adic integer s (33 = 1 mod 8), so L has two places over 2. At each, the
# completion is Q_2(sqrt3): ramified over Q_2, residue field F_2, uniformiser 1 + sqrt3, and sqrt3 is a
# unit with sqrt3 = 1 mod (1 + sqrt3). Neither -1 nor -3 is a square in Q_2, so i is not in Q_2(sqrt3),
# and K_w = Q_2(sqrt3)(i) = Q_2(sqrt3)(sqrt-3) is the unramified quadratic extension: the place is
# INERT in K/L, with O_w = Z_2[sqrt3][w] (w = (-1 + sqrt-3)/2) and residue field F_4.
# A unit vector u (u ubar = 1) therefore has |u|_w = 1, so u is in O_w^x and its residue is one of
# the three nonzero elements of F_4. Colour z by the residue of z - rep(z), where rep(z) is the
# 2-adic fractional part of z's coordinates on the Z_2-basis 1, sqrt3, w, sqrt3 w of O_w: adding u
# adds its (nonzero) residue, so the colouring is proper with the 4 colours of F_4.
# (Moorhouse 2010 left this field open; Madore, arXiv 1509.07023, proved 4 <= chi <= 5.)

def _two_adic_scaled(r0, r1, s, K, N):
    """(r0 + r1 s) * 2^K mod 2^N as an integer, for rationals r0, r1 whose 2-adic valuation is >= -K."""
    tot = 0
    for r, m in ((Fr(r0), 1), (Fr(r1), s)):
        if r == 0:
            continue
        n, d = r.numerator, r.denominator
        k = 0
        while d % 2 == 0:
            d //= 2
            k += 1
        assert k <= K, "denominator too 2-divisible for the chosen scale"
        tot += n * m * (1 << (K - k)) * pow(d, -1, 1 << N)
    return tot % (1 << N)


def q311_coordinates(p, place=1, K=64):
    """z = x + iy on the Z_2-basis 1, sqrt3, w, sqrt3 w of O_w, at the place sqrt33 -> place * s.

    Returns four integers c_j = coefficient_j * 2^K mod 2^PRECISION. The coefficient is a 2-adic
    integer iff c_j = 0 mod 2^K, and its units digit is bit K of c_j."""
    F = p.x.field
    X = dict(zip(F._prod, p.x.c)); Y = dict(zip(F._prod, p.y.c))
    assert all(q in (1, 3, 11, 33) for q, v in list(X.items()) + list(Y.items()) if v), "not in Q(sqrt3, sqrt11)"
    s = place * _S
    # sqrt11 = sqrt33 / sqrt3 = (s / 3) sqrt3, so x = xa + xb sqrt3 with xa = x1 + x33 s, xb = x3 + x11 s / 3
    xa = (Fr(X.get(1, 0)), Fr(X.get(33, 0))); xb = (Fr(X.get(3, 0)), Fr(X.get(11, 0)) / 3)
    ya = (Fr(Y.get(1, 0)), Fr(Y.get(33, 0))); yb = (Fr(Y.get(3, 0)), Fr(Y.get(11, 0)) / 3)
    # i = sqrt-3 / sqrt3 = (1 + 2w) sqrt3 / 3, so z = (xa + yb) + (xb + ya/3) sqrt3 + (2 yb) w + (2 ya/3) sqrt3 w
    coeffs = [(xa[0] + yb[0], xa[1] + yb[1]), (xb[0] + ya[0] / 3, xb[1] + ya[1] / 3),
              (2 * yb[0], 2 * yb[1]), (2 * ya[0] / 3, 2 * ya[1] / 3)]
    return [_two_adic_scaled(r0, r1, s, K, PRECISION) for r0, r1 in coeffs]


def q311_colour(p, place=1, K=64):
    """The colour 0..3 of a point of Q(sqrt3, sqrt11)^2: the residue in F_4 = F_2[w] of z - rep(z).

    Coordinates are on the basis 1, sqrt3, sqrt11, sqrt33 (hn.field.Field((3, 11)) or a subfield of
    a larger Field). sqrt3 = 1 mod (1 + sqrt3), so the residue of a + b sqrt3 + c w + d sqrt3 w is
    (a + b) + (c + d) w with a, b, c, d the units digits."""
    a, b, c, d = ((v >> K) & 1 for v in q311_coordinates(p, place, K))
    return ((a + b) & 1) + 2 * ((c + d) & 1)


# ---------------------------------------------------------------------------------------------
# The plane Q(sqrt2, sqrt3)^2 is 4-colourable too, so its chromatic number is 4.
#
# L = Q(sqrt2, sqrt3) has a single place v over 2: in Q_2* / squares, 2, 3 = -5 and 6 = -10 generate a
# group of order 4, so L_v = Q_2(sqrt2, sqrt3) has degree 4. That group misses 5 = -3, the class of
# the unramified quadratic extension, so v is totally ramified with residue field F_2, and
# v(c) = v_2(N_{L/Q}(c)). It also misses -1, so i is not in L_v and v does not split in K = L(i).
# K_w contains sqrt-3 = i sqrt3, so K_w = L_v(w) with w = (-1 + sqrt-3)/2: unramified over L_v, with
# residue field F_4 and O_w = O_v + O_v w. With i = (2w + 1)/sqrt3, z = x + iy = a + b w where
# a = x + y/sqrt3 and b = 2y/sqrt3. A unit vector u (u ubar = 1) has |u|_w = 1, so its a and b are
# integral and not both in the maximal ideal: its residue in F_4 = F_2[w] is nonzero. Colouring z by
# the residue of z - rep(z) is therefore proper. Voronov (Polymath16, 2021) conjectured
# chi(Q(i, sqrt2, sqrt3)) = 4 together with the case Q(i, sqrt3, sqrt11) above.

def q23_norm(c):
    """N_{L/Q}(c) for c in Q(sqrt2, sqrt3), coefficients on the basis 1, sqrt2, sqrt3, sqrt6."""
    F = c.field
    p, q, r, s = c.c
    conj = [F.element([p, e2 * q, e3 * r, e2 * e3 * s]) for e2, e3 in ((-1, 1), (1, -1), (-1, -1))]
    n = c * conj[0] * conj[1] * conj[2]
    assert all(x == 0 for x in n.c[1:])
    return n.c[0]


def q23_valuation(c):
    """The valuation of c at the place of Q(sqrt2, sqrt3) over 2, normalised so a uniformiser has 1."""
    q = Fr(q23_norm(c))
    if q == 0:
        return None
    k, n, d = 0, q.numerator, q.denominator
    while n % 2 == 0:
        n //= 2
        k += 1
    while d % 2 == 0:
        d //= 2
        k -= 1
    return k


def q23_ab(p):
    """z = x + iy written as a + b w, w = (-1 + sqrt-3)/2: a = x + y sqrt3/3, b = 2 y sqrt3/3."""
    F = p.x.field
    s3 = F.sqrt(3)
    return p.x + p.y * s3 * F.rational(Fr(1, 3)), p.y * s3 * F.rational(Fr(2, 3))


def q23_colour(p, base):
    """The colour 0..3 of p relative to a base point of its connected component: the residue in F_4 of
    z - z_base, whose a and b are integral when p and base are joined by unit steps."""
    a, b = q23_ab(p)
    a0, b0 = q23_ab(base)
    va, vb = q23_valuation(a - a0), q23_valuation(b - b0)
    assert (va is None or va >= 0) and (vb is None or vb >= 0), "not in the same coset of O_w"
    return (1 if va == 0 else 0) + 2 * (1 if vb == 0 else 0)
