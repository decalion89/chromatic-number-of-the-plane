#!/usr/bin/env python3
"""Independent referee check of Lemma 13 (lem:seven) of three-colours.tex.

Written from scratch; uses only the Python standard library.  Everything is exact:
F_49 = Z[i]/7 as pairs (a, b) mod 7, complex numbers as pairs of Fractions.
Each check prints PASS/FAIL; the script exits with status 1 if any check fails.
"""
from fractions import Fraction as Fr
from itertools import product
import sys

FAILS = []


def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""))
    if not cond:
        FAILS.append(name)


P = 7
# ---------------------------------------------------------------- F_49 = Z[i]/7
F49 = [(a, b) for a in range(P) for b in range(P)]
ZERO = (0, 0)


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def mul(x, y):
    return ((x[0] * y[0] - x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)


def conj(x):
    return (x[0] % P, (-x[1]) % P)


def norm(x):  # x * conj(x), an element of F_7
    return (x[0] * x[0] + x[1] * x[1]) % P


def re(x):
    return x[0] % P


def power(x, n):
    r = (1, 0)
    for _ in range(n):
        r = mul(r, x)
    return r


SQUARES = {(s * s) % P for s in range(1, P)}
check("nonzero squares mod 7 are {1,2,4}", SQUARES == {1, 2, 4}, str(sorted(SQUARES)))
check("-1 is not a square mod 7", (P - 1) not in SQUARES)

# mu_8 = elements of norm one
MU8 = [z for z in F49 if norm(z) == 1]
check("|mu_8| = 8", len(MU8) == 8, str(MU8))
check("mu_8 = roots of z^8 - 1 (norm one <=> z^8 = 1)",
      set(MU8) == {z for z in F49 if power(z, 8) == (1, 0)})
check("conjugation is Frobenius z -> z^7 on F_49", all(conj(z) == power(z, 7) for z in F49))
rho_t = (2, 5)  # rho~ = 2+5i
check("rho = (3+4i)/5 = 2+5i mod 7  (5*(2+5i) = 3+4i mod 7)", mul((5, 0), rho_t) == (3, 4))
check("rho~ has norm 1 (29 = 1 mod 7)", norm(rho_t) == 1 and (2 * 2 + 5 * 5) == 29)
check("rho~^2 = -21+20i (as integers)", (2 * 2 - 5 * 5, 2 * 2 * 5) == (-21, 20))
check("rho~^2 = -i mod 7", mul(rho_t, rho_t) == (0, P - 1))
order = next(n for n in range(1, 49) if power(rho_t, n) == (1, 0))
check("rho~ has order 8 (generates mu_8)", order == 8 and {power(rho_t, j) for j in range(8)} == set(MU8))
check("conj(rho~) = rho~^{-1} mod 7", mul(rho_t, conj(rho_t)) == (1, 0))

# ------------------------------------------------------------------ (i)
A7_def = {e for e in F49 if all(re(mul(conj(e), z)) in {2, 3, 4, 5} for z in MU8)}
LIST = {(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)}
check("(i) A'_7 (from its definition) equals the listed 8 elements", A7_def == LIST, str(sorted(A7_def)))
check("(i) A'_7 = {e : e ebar = -1}", A7_def == {e for e in F49 if norm(e) == P - 1})
check("(i) Re(ebar z) = ax + by",
      all(re(mul(conj(e), z)) == (e[0] * z[0] + e[1] * z[1]) % P for e in F49 for z in F49))
check("(i) solutions of a^2+b^2=-1 have {a^2,b^2}={2,4}",
      all({(a * a) % P, (b * b) % P} == {2, 4} for (a, b) in A7_def))

# line / conic argument, exhaustively for every e != 0 and every t
ok_param, ok_formula, ok_values = True, True, True
for e in F49:
    if e == ZERO:
        continue
    a, b = e
    N = norm(e)
    Ninv = pow(N, P - 2, P)
    values = {re(mul(conj(e), z)) for z in MU8}
    for t in range(P):
        line = {(x, y) for x in range(P) for y in range(P) if (a * x + b * y) % P == t}
        param = {((t * Ninv * a - s * b) % P, (t * Ninv * b + s * a) % P) for s in range(P)}
        if line != param or len(param) != 7:
            ok_param = False
        for s in range(P):
            pt = ((t * Ninv * a - s * b) % P, (t * Ninv * b + s * a) % P)
            on_conic = norm(pt) == 1
            if on_conic != ((N * N * s * s) % P == (N - t * t) % P):
                ok_formula = False
        is_value = t in values
        if is_value != (((N - t * t) % P == 0) or ((N - t * t) % P in SQUARES)):
            ok_values = False
check("(i) the line ax+by=t is {t/(e ebar)(a,b) + s(-b,a)} (7 points), all e!=0, t", ok_param)
check("(i) such a point is on x^2+y^2=1 iff (e ebar)^2 s^2 = e ebar - t^2, all e!=0, t, s", ok_formula)
check("(i) t is a value of Re(ebar z) on mu_8 iff e ebar - t^2 is 0 or a square, all e!=0, t", ok_values)
check("(i) e in A'_7 iff e!=0 and e ebar, e ebar - 1 both in {3,5,6}",
      A7_def == {e for e in F49 if e != ZERO and norm(e) in {3, 5, 6} and (norm(e) - 1) % P in {3, 5, 6}})

e0 = (2, 3)
check("(i) A'_7 = e0*mu_8 for every e0 in A'_7 (single mu_8-orbit)",
      all({mul(x, z) for z in MU8} == A7_def for x in A7_def))
check("(i) A'_7 is closed under conjugation", {conj(e) for e in A7_def} == A7_def)
s = ZERO
for e in A7_def:
    s = add(s, e)
check("(i) the elements of A'_7 add up to 0", s == ZERO)
smu = ZERO
for z in MU8:
    smu = add(smu, z)
check("(i) sum of mu_8 is 0", smu == ZERO)
check("(i) 0 not in A'_7", ZERO not in A7_def)

# ------------------------------------------------------------------ (ii)
# all additive subgroups of F_7^2: {0}, the 8 lines through 0, the whole group
subgroups = set()
for d in F49:
    subgroups.add(frozenset(((s * d[0]) % P, (s * d[1]) % P) for s in range(P)))
subgroups.add(frozenset(F49))
nonzero_subgroups = [H for H in subgroups if len(H) > 1]
check("(ii) number of nonzero additive subgroups is 9 (8 lines + whole)", len(nonzero_subgroups) == 9)
bad = [(H, c) for H in nonzero_subgroups for c in F49 if {add(c, h) for h in H} <= A7_def]
check("(ii) no coset of a nonzero additive subgroup lies in A'_7", not bad)
max_on_line = max(len({add(c, h) for h in H} & A7_def) for H in nonzero_subgroups if len(H) == 7 for c in F49)
check("(ii) every affine line meets A'_7 in at most 2 points", max_on_line <= 2, "max = %d" % max_on_line)
# the degree-2 polynomial in s, as in the proof
ok_poly = True
for e in F49:
    for d in F49:
        if d == ZERO:
            continue
        tr = (2 * (e[0] * d[0] + e[1] * d[1])) % P  # e dbar + ebar d = 2 Re(e dbar)
        for s_ in range(P):
            lhs = norm(add(e, mul((s_, 0), d)))
            rhs = (norm(e) + s_ * tr + s_ * s_ * norm(d)) % P
            if lhs != rhs:
                ok_poly = False
check("(ii) N(e+sd) = e ebar + s(e dbar + ebar d) + s^2 d dbar, all e, d!=0, s", ok_poly)

# ------------------------------------------------------------------ (iii)
check("(iii) Re(A'_7) = {2,3,4,5}", {re(e) for e in A7_def} == {2, 3, 4, 5})
check("both coordinates of every element of A'_7 lie in {2,3,4,5}",
      all(a in {2, 3, 4, 5} and b in {2, 3, 4, 5} for (a, b) in A7_def))
check("each value 2,3,4,5 is attained by Re(ebar z) on mu_8, for every e in A'_7 (least margin 2/7)",
      all({re(mul(conj(e), z)) for z in MU8} == {2, 3, 4, 5} for e in A7_def))

# ------------------------------------------------------------------ (iv): exact complex arithmetic
class C:
    __slots__ = ("x", "y")

    def __init__(self, x, y=0):
        self.x, self.y = Fr(x), Fr(y)

    def __add__(s, o):
        return C(s.x + o.x, s.y + o.y)

    def __sub__(s, o):
        return C(s.x - o.x, s.y - o.y)

    def __mul__(s, o):
        return C(s.x * o.x - s.y * o.y, s.x * o.y + s.y * o.x)

    def __eq__(s, o):
        return s.x == o.x and s.y == o.y

    def __hash__(s):
        return hash((s.x, s.y))

    def conj(s):
        return C(s.x, -s.y)

    def frac_centered(s):
        """representative of s mod Z[i] with coordinates in [-1/2, 1/2)"""
        def fc(t):
            r = t - (t.numerator // t.denominator)
            return r - 1 if r >= Fr(1, 2) else r
        return C(fc(s.x), fc(s.y))

    def __repr__(s):
        if not s.y:
            return "%s" % s.x
        return "(%s %s %s i)" % (s.x, "+" if s.y > 0 else "-", abs(s.y))


I = C(0, 1)
MINUS_I = C(0, -1)
RHO = C(Fr(3, 5), Fr(4, 5))
H = C(Fr(1, 2), Fr(1, 2))
RHO_T = C(2, 5)
ALLOWED = {Fr(1, 14), Fr(-1, 14), Fr(3, 14), Fr(-3, 14)}


def gpow(z, n):
    r = C(1)
    for _ in range(n):
        r = r * z
    return r


def zeta(e, j):
    """zeta_j(eps) for 7 eps = e (a Gaussian integer); rho~^j read mod 8"""
    eps = C(Fr(e[0], 7), Fr(e[1], 7))
    w = eps.conj() * gpow(RHO_T, j % 8)
    z = (w - H).frac_centered()
    assert z.x in ALLOWED and z.y in ALLOWED, ("zeta undefined", e, j, z)
    return z


def nu(e, j):
    return zeta(e, j) - RHO * zeta(e, j - 1)


JR = range(-24, 25)
ok_def = True
try:
    for e in A7_def:
        for j in JR:
            zeta(e, j)
except AssertionError:
    ok_def = False
check("(iv) zeta_j(eps) is defined (in {+-1/14,+-3/14}^2) for all 8 classes and -24<=j<=24", ok_def)
check("(iv) zeta_{j+2} = -i zeta_j for all classes, -24<=j<=22",
      all(zeta(e, j + 2) == MINUS_I * zeta(e, j) for e in A7_def for j in range(-24, 23)))
check("(iv) {+-1/14,+-3/14}^2 is stable under multiplication by -i",
      all((MINUS_I * C(x, y)).x in ALLOWED and (MINUS_I * C(x, y)).y in ALLOWED for x in ALLOWED for y in ALLOWED))
check("(iv) -i h = h mod Z[i]", (MINUS_I * H - H).frac_centered() == C(0))
check("(iv) nu_{j+2} = -i nu_j for all classes", all(nu(e, j + 2) == MINUS_I * nu(e, j) for e in A7_def for j in range(-22, 23)))

# hand computations for 7 eps = 2+3i
check("(iv) conj(2+3i) = 2+4i mod 7", conj((2, 3)) == (2, 4))
check("(iv) (2+4i) rho~ = 5+4i mod 7", mul((2, 4), rho_t) == (5, 4))
check("(iv) (2+4i) rho~^2 = 4+5i mod 7", mul((2, 4), mul(rho_t, rho_t)) == (4, 5))
check("(iv) (2+4i)(2+5i) = -16+18i exactly", C(2, 4) * C(2, 5) == C(-16, 18))
z0, z1, z2 = zeta((2, 3), 0), zeta((2, 3), 1), zeta((2, 3), 2)
check("(iv) zeta_0 = (-3+i)/14", z0 == C(Fr(-3, 14), Fr(1, 14)), repr(z0))
check("(iv) zeta_1 = (3+i)/14", z1 == C(Fr(3, 14), Fr(1, 14)), repr(z1))
check("(iv) zeta_2 = (1+3i)/14", z2 == C(Fr(1, 14), Fr(3, 14)), repr(z2))
check("(iv) rho zeta_0 = -(13+9i)/70", RHO * z0 == C(Fr(-13, 70), Fr(-9, 70)), repr(RHO * z0))
check("(iv) rho zeta_1 = (1+3i)/14", RHO * z1 == C(Fr(1, 14), Fr(3, 14)), repr(RHO * z1))
check("(iv) nu_1 = (2+i)/5", nu((2, 3), 1) == C(Fr(2, 5), Fr(1, 5)))
check("(iv) nu_2 = 0", nu((2, 3), 2) == C(0))
BASE = C(Fr(2, 5), Fr(1, 5))
ok_pat = all(nu((2, 3), 2 * m) == C(0) and nu((2, 3), 2 * m + 1) == gpow(MINUS_I, m % 4) * BASE
             for m in range(-12, 12))
check("(iv) nu_{2m} = 0 and nu_{2m+1} = (-i)^m (2+i)/5, -12<=m<12", ok_pat)

UNITS = [C(1), I, C(-1), MINUS_I]
DIG7 = {u * BASE for u in UNITS}
ok_form = True
for e in A7_def:
    seq = [nu(e, j) for j in range(-16, 17)]
    if any(seq[k] != seq[k + 8] for k in range(len(seq) - 8)):
        ok_form = False
    zero_pos = [k for k in range(len(seq)) if seq[k] == C(0)]
    if not all(((k - zero_pos[0]) % 2 == 0) for k in zero_pos) or len(zero_pos) * 2 not in (len(seq) - 1, len(seq) + 1):
        ok_form = False
    if any(x != C(0) and x not in DIG7 for x in seq):
        ok_form = False
check("(iv) every class: period 8, 0 alternates with digits u(2+i)/5, u in mu_4", ok_form)
check("(iv) exact period is 8 (not 4) for every class",
      all(any(nu(e, j) != nu(e, j + 4) for j in range(8)) for e in A7_def))

# shift rule
ok_shift = all(zeta(mul(e, power(conj(rho_t), j % 8)), t) == zeta(e, t + j)
               for e in A7_def for j in range(-8, 9) for t in range(-10, 11))
check("(iv) zeta_t(eps') = zeta_{t+j}(eps) when 7eps' = 7eps conj(rho~)^j mod 7", ok_shift)
shift_of = {}
for e in A7_def:
    js = [j for j in range(8) if mul((2, 3), power(conj(rho_t), j)) == e]
    shift_of[e] = js
check("(iv) every class is (2+3i)conj(rho~)^j for exactly one j mod 8 (8 classes <-> 8 shifts)",
      all(len(v) == 1 for v in shift_of.values()) and sorted(v[0] for v in shift_of.values()) == list(range(8)))
seqs7 = {e: tuple(nu(e, j) for j in range(0, 8)) for e in A7_def}
check("(iv) the 8 classes give 8 distinct digit sequences", len(set(seqs7.values())) == 8)

# ------------------------------------------------------------- type-point digits from strip indices
def floor_fr(t):
    return t.numerator // t.denominator


def rho_pow(t):
    return gpow(RHO, t) if t >= 0 else gpow(RHO.conj(), -t)


def digits_of_point(c, ts):
    """nu_t(c) = g_t - rho g_{t-1}, g_t = cbar rho^t - h - n_t, n_t = floors of the two coordinates
    (n_{rho^t} = floor Re(cbar rho^t), n_{-i rho^t} = floor Re(cbar (-i) rho^t) = floor Im(cbar rho^t))."""
    def g(t):
        w = c.conj() * rho_pow(t)
        n1 = floor_fr((w).x)
        n2 = floor_fr((c.conj() * MINUS_I * rho_pow(t)).x)
        check_ok = (n2 == floor_fr(w.y))
        assert check_ok
        return w - H - C(n1, n2)
    return tuple(g(t) - RHO * g(t - 1) for t in ts)


def q_seq(eta0, ts):
    # digits of the proof of Lemma 12: nu_t = -(1+3i)/10 * eta_{t-1}, eta_{t+1} = -i eta_t
    k = C(Fr(-1, 10), Fr(-3, 10))
    return tuple(k * gpow(MINUS_I, (t - 1) % 4) * eta0 for t in ts)


ETAS = [C(1, 1), C(1, -1), C(-1, 1), C(-1, -1)]
TS = (-1, 0, 1, 2)
# level 325 (Proposition 9 / Lemma 14)
ok_c = digits_of_point(C(Fr(325, 2), Fr(325, 2)), TS) == (C(0),) * 4
check("type point 325h has digits 0,0,0,0 at t=-1..2", ok_c)
qdig = []
okq = True
for al in (1, 2):
    for be in (1, 2):
        d = digits_of_point(C(Fr(325 * al, 3), Fr(325 * be, 3)), TS)
        nu0 = d[0]
        if not (nu0 != C(0) and d == (nu0, MINUS_I * nu0, C(-1) * nu0, I * nu0)):
            okq = False
        if d not in {q_seq(eta, TS) for eta in ETAS}:
            okq = False
        qdig.append(d)
check("type-q points (325/3)(a+bi): digits (nu,-i nu,-nu,i nu), nu!=0, = digits of the proof of Lemma 12",
      okq and len(set(qdig)) == 4)
ok7 = True
n7 = 0
for a, b in product(range(7), range(7)):
    e = mul((3, 0), (a, b))
    if e not in A7_def:
        continue
    n7 += 1
    d = digits_of_point(C(Fr(325 * a, 7), Fr(325 * b, 7)), TS)
    if d != tuple(nu(e, t) for t in TS):
        ok7 = False
check("type-7 points (325/7)(a+bi), 3(a+bi) in A'_7: 8 classes, digits = zeta-digits of eps, 7eps=3(a+bi)",
      ok7 and n7 == 8)
check("325 = 3 mod 7 and 325 = 1 mod 3", 325 % 7 == 3 and 325 % 3 == 1)
# level 65 (Proposition 8 / Lemma 12) for comparison
ok65 = digits_of_point(C(Fr(65, 2), Fr(65, 2)), (0, 1)) == (C(0), C(0))
for al in (1, 2):
    for be in (1, 2):
        d = digits_of_point(C(Fr(65 * al, 3), Fr(65 * be, 3)), (0, 1))
        if not (d[0] != C(0) and d[1] == MINUS_I * d[0] and d in {q_seq(eta, (0, 1)) for eta in ETAS}):
            ok65 = False
check("level 65: type c digits (0,0); type q digits (nu,-i nu) of the q family", ok65)
pts1046 = [(1, 2), (1, 5), (2, 1), (2, 6), (5, 1), (5, 6), (6, 2), (6, 5)]
check("the 7-torsion points (65/7)(a+bi) listed after Prop. 8 are exactly those with 2(a+bi) in A'_7",
      {(a, b) for a in range(7) for b in range(7) if mul((2, 0), (a, b)) in A7_def} == set(pts1046))
check("lambda(a+bi)=2a+3b corresponds to (65/7)(1+5i): 2(1+5i) = 2+3i mod 7", mul((2, 0), (1, 5)) == (2, 3))

# ------------------------------------------------------------- two consecutive digits determine type & class
FAMILY = {("c", None): lambda t: C(0)}
for eta in ETAS:
    FAMILY[("q", (eta.x, eta.y))] = (lambda eta_: (lambda t: q_seq(eta_, (t,))[0]))(eta)
for e in A7_def:
    FAMILY[("7", e)] = (lambda e_: (lambda t: nu(e_, t)))(e)
ok_inj = True
for t in range(-12, 13):
    seen = {}
    for key, f in FAMILY.items():
        pair = (f(t), f(t + 1))
        if pair in seen:
            ok_inj = False
        seen[pair] = key
check("(iv) for every position t in [-12,12], (nu_t, nu_{t+1}) determines type and class among the 13 sequences",
      ok_inj and len(FAMILY) == 13)
# window chaining in Lemma 14: window at j+1 of class eps_{j+1}, 7eps_{j+1} = 7eps_j conj(rho~)
ok_chain = all(tuple(nu(mul(e, conj(rho_t)), t) for t in TS) == tuple(nu(e, t + 1) for t in TS) for e in A7_def)
check("Lemma 14: the type-7 window at j+1 has class 7eps conj(rho~) (its digits are the shifted ones)", ok_chain)

# ------------------------------------------------------------- facts used in Lemma 14
ok_N = True
for k in (6, 12):
    N = 5 ** k
    check("5^%d = 1 mod 7" % k, N % 7 == 1)
    for j in range(-k, k + 1):
        w = C(N) * rho_pow(j)
        if w.x.denominator != 1 or w.y.denominator != 1:
            ok_N = False
            continue
        r = gpow(RHO_T, j % 8)
        if (int(w.x) - int(r.x)) % 7 or (int(w.y) - int(r.y)) % 7:
            ok_N = False
check("N rho^j in Z[i] and N rho^j = rho~^j mod 7 for N=5^6, 5^12, |j|<=k", ok_N)
rprime = (Fr(1, 4) + Fr(3, 14)) ** 2 * 2  # r'^2
check("r' = (1/4+3/14) sqrt2 < 0.66", rprime < Fr(66, 100) ** 2, "r'^2 = %s" % rprime)
check("|zeta| <= 3 sqrt2/14 for all zeta values",
      all(z.x ** 2 + z.y ** 2 <= 2 * Fr(3, 14) ** 2 for e in A7_def for z in [zeta(e, j) for j in range(8)]))

# ------------------------------------------------------------- the step of Lemma 17 (B3seven), exhaustively
bad17 = [(c, w) for c in F49 for w in F49 if w != ZERO and {add(c, mul(o, w)) for o in MU8} <= A7_def and c != ZERO]
check("Lemma 17 core: c + mu_8 w in A'_7 with w != 0 forces c = 0 (all c, w)", not bad17)

print()
print("FAILED:", FAILS if FAILS else "none")
sys.exit(1 if FAILS else 0)
