"""Referee's own exact checks for Section 10 (places above 7, type points, digit words).

Everything is recomputed from the definitions in the text, with Fractions / integers only.
"""
from fractions import Fraction as Fr
from itertools import product

P = 7


# ---------------------------------------------------------------- F_49 = Z[i]/7
def fmul(x, y):
    return ((x[0] * y[0] - x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)


def fconj(x):
    return (x[0] % P, (-x[1]) % P)


F49 = [(a, b) for a in range(P) for b in range(P)]
mu8 = [z for z in F49 if fmul(z, fconj(z)) == (1, 0)]
assert len(mu8) == 8
# A'_7 from the definition: Re(conj(e) z) in {2,3,4,5} for every z in mu_8
A7 = sorted(e for e in F49 if all(fmul(fconj(e), z)[0] in (2, 3, 4, 5) for z in mu8))
LIST = sorted([(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)])
print("Lemma 13(i) list:", A7 == LIST, A7)
orbit = sorted({fmul(A7[0], z) for z in mu8})
print("  single mu_8-orbit:", orbit == A7)
print("  closed under conjugation:", sorted(fconj(e) for e in A7) == A7)
print("  sum zero:", (sum(e[0] for e in A7) % P, sum(e[1] for e in A7) % P) == (0, 0), " 0 not in A'_7:", (0, 0) not in A7)
# (ii) cosets of nonzero additive subgroups: the subgroups of F_7^2 are 0, the 8 lines, and everything
lines = []
for d in [(1, t) for t in range(P)] + [(0, 1)]:
    lines.append({((t * d[0]) % P, (t * d[1]) % P) for t in range(P)})
bad = []
for L in lines:
    for x in F49:
        coset = {((x[0] + y[0]) % P, (x[1] + y[1]) % P) for y in L}
        if coset <= set(A7):
            bad.append((x, L))
print("Lemma 13(ii) no coset of a nonzero subgroup (49 > 8 handles the whole group):", not bad)
print("Lemma 13(iii) Re(A'_7) =", sorted({e[0] for e in A7}))

# ---------------------------------------------------------------- complex numbers with Fraction coordinates
def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cdiv(x, y):
    n = y[0] * y[0] + y[1] * y[1]
    return (Fr(x[0] * y[0] + x[1] * y[1], n), Fr(x[1] * y[0] - x[0] * y[1], n))


def cpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z, e = cdiv((1, 0), z), -e
    for _ in range(e):
        w = cmul(w, z)
    return w


def conj(x):
    return (x[0], -x[1])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def scal(c, x):
    return (c * x[0], c * x[1])


def fl(q):
    return q.numerator // q.denominator


RHO = cdiv((2, 1), (2, -1))
SIG = cdiv((3, 2), (3, -2))
assert RHO == (Fr(3, 5), Fr(4, 5)) and SIG == (Fr(5, 13), Fr(12, 13))
H = (Fr(1, 2), Fr(1, 2))
I = (Fr(0), Fr(1))

# rho~ = 2+5i = rho mod 7, order 8, rho~^2 = -i
RT = (2, 5)
print("rho = 2+5i mod 7:", ((3 * 3) % 7, (3 * 4) % 7) == RT, " (5^-1 = 3 mod 7)")
pw = [(1, 0)]
for _ in range(8):
    pw.append(fmul(pw[-1], RT))
print("rho~ has order 8, rho~^2 = -i:", pw[8] == (1, 0) and all(pw[t] != (1, 0) for t in range(1, 8)), pw[2] == (0, 6))


def rt_pow(j):
    return pw[j % 8]


def frac_part(x):
    return (x[0] - fl(x[0]), x[1] - fl(x[1]))


def g_of(w):
    """w = h + g + n with n = floors (strip indices of Re and Im), g in [-1/2,1/2)^2"""
    n = (fl(w[0]), fl(w[1]))
    return (w[0] - H[0] - n[0], w[1] - H[1] - n[1]), n


# zeta_j(eps) for eps = e/7, e in A'_7: conj(eps) rho~^j in h + zeta + Z[i]
ZSET = {Fr(1, 14), Fr(-1, 14), Fr(3, 14), Fr(-3, 14)}


def zeta(e, j):
    ce = fconj(e)
    w = fmul(ce, rt_pow(j))
    w = (Fr(w[0], 7), Fr(w[1], 7))
    g, _ = g_of(w)
    assert g[0] in ZSET and g[1] in ZSET, (e, j, g)
    return g


def digit7(e, j):
    return sub(zeta(e, j), cmul(RHO, zeta(e, j - 1)))


seq7 = {e: [digit7(e, j) for j in range(16)] for e in A7}
ok_period = all(seq7[e][j] == seq7[e][j + 8] for e in A7 for j in range(8))
U4 = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))]
D7 = [cmul(u, (Fr(2, 5), Fr(1, 5))) for u in U4]
ZERO = (Fr(0), Fr(0))
alt = True
for e in A7:
    s = seq7[e][:8]
    par = [d == ZERO for d in s]
    if not (all(par[j] != par[j + 1] for j in range(7))):
        alt = False
    if not all(d == ZERO or d in D7 for d in s):
        alt = False
    if sorted(d for d in s if d != ZERO) != sorted(D7):
        alt = False
print("Lemma 13(iv): period 8:", ok_period, "; 0 alternates with u(2+i)/5, each u once per period:", alt)
S0 = seq7[A7[0]][:8]
shifts = {e: [t for t in range(8) if all(seq7[e][j] == S0[(j + t) % 8] for j in range(8))] for e in A7}
print("  the 8 classes are the 8 shifts of one sequence:", sorted(sum(shifts.values(), [])) == list(range(8)), shifts)
# relation between classes: e*rho~^j shifts by j
print("  e*rho~ shifts the sequence by one:",
      all(all(seq7[fmul(e, RT)][j + 1] == seq7[e][j] for j in range(8)) for e in A7))

# type q digits (Lemma 12): nu = -(1+3i) eta / 10, nu_{t+1} = -i nu_t
ETAS = [(Fr(a), Fr(b)) for a in (1, -1) for b in (1, -1)]
Dq = [scal(Fr(-1, 10), cmul((Fr(1), Fr(3)), eta)) for eta in ETAS]
print("  type-q digits = type-7 nonzero digits as sets:", sorted(Dq) == sorted(D7))
seqq = {}
for nu in Dq:
    s = [nu]
    for _ in range(15):
        s.append(cmul(scal(-1, I), s[-1]))
    seqq[nu] = s
seqc = [ZERO] * 16

# all (type, shift) sequences, as functions of position (period 8 suffices: q has period 4)
SEQS = [("c", 0, seqc)] + [("q", k, seqq[nu]) for k, nu in enumerate(sorted(Dq))] + \
       [("7", t, [S0[(j + t) % 8] for j in range(16)]) for t in range(8)]


def pairs_determine(L):
    seen = {}
    for name, sh, s in SEQS:
        for p in range(8):
            key = tuple(s[p:p + L])
            # the "sequence aligned at p" = (name, shifted sequence)
            aligned = (name, tuple(s[(p + j) % 8] for j in range(8)) if name != "c" else ())
            if key in seen and seen[key] != aligned:
                return False
            seen[key] = aligned
    return True


print("  two consecutive digits determine type and aligned sequence:", pairs_determine(2),
      "; one digit:", pairs_determine(1))

# ------------------------------------------------- overlap automaton for Lemma 14 (windows of 4, overlap 3)
W4 = set()
for name, sh, s in SEQS:
    for p in range(8):
        W4.add(tuple(s[p:p + 4]))
print("number of allowed 4-digit windows (c, q, 7):", len(W4))


def pure(word):
    for name, sh, s in SEQS:
        for p in range(8):
            if all(word[j] == s[(p + j) % 8] for j in range(len(word))):
                return True
    return False


# all words of length Lw whose every 4-window is allowed
ok = True
for Lw in range(4, 15):
    words = [w for w in W4]
    for _ in range(Lw - 4):
        words = [w + (d,) for w in words for d in set(x[3] for x in W4) if w[-3:] + (d,) in W4]
    if not all(pure(w) for w in words):
        ok = False
print("every word (length 4..14) all of whose 4-windows are allowed is a segment of one sequence:", ok)

# ------------------------------------------------- type points at modulus 65 and 325: digits
def digits_at(c, ts):
    """nu_t(c) for t in ts, from floors of Re/Im of conj(c) rho^t"""
    res = {}
    for t in ts:
        g1, _ = g_of(cmul(conj(c), cpow(RHO, t)))
        g0, _ = g_of(cmul(conj(c), cpow(RHO, t - 1)))
        res[t] = sub(g1, cmul(RHO, g0))
    return res


def typepoints(N):
    pts = [("c", None, scal(Fr(N), H))]
    for a, b in product((1, 2), repeat=2):
        pts.append(("q", (a, b), (Fr(N * a, 3), Fr(N * b, 3))))
    return pts


# modulus 65: (nu_0, nu_1) at the type points
for name, ab, c in typepoints(65):
    d = digits_at(c, [0, 1])
    if name == "c":
        assert d[0] == ZERO and d[1] == ZERO
    else:
        assert d[0] in Dq and d[1] == cmul(scal(-1, I), d[0]), (ab, d)
print("modulus 65: type points give (0,0) or (nu,-i nu) with nu in the q-digits: True")

# modulus 325: (nu_-1..nu_2)
for name, ab, c in typepoints(325):
    d = digits_at(c, [-1, 0, 1, 2])
    w = [d[t] for t in (-1, 0, 1, 2)]
    if name == "c":
        assert all(x == ZERO for x in w)
    else:
        nu = w[0]
        mi = scal(-1, I)
        assert nu in Dq and w[1] == cmul(mi, nu) and w[2] == scal(-1, nu) and w[3] == cmul(I, nu), (ab, w)
print("modulus 325: type c gives (0,0,0,0), type q gives (nu,-i nu,-nu,i nu): True")

# type 7 points (325/7)(a+bi): compare three descriptions
G21 = [cmul(cpow(RHO, j), cpow(SIG, l)) for j in range(-2, 3) for l in range(-1, 2)]


def all_values_in_2to5(c, G):
    for g in G:
        w = cmul(conj(c), g)
        for x in (w[0], w[1]):
            f = x - fl(x)
            if f not in (Fr(2, 7), Fr(3, 7), Fr(4, 7), Fr(5, 7)):
                return False
    return True


sel_text, sel_vals, sel_zeta = set(), set(), set()
for a, b in product(range(7), repeat=2):
    c = (Fr(325 * a, 7), Fr(325 * b, 7))
    e = ((3 * a) % 7, (3 * b) % 7)
    if e in A7:
        sel_text.add((a, b))
        # values at rho^t, |t| <= 2, are those of eps = e/7: conj(c) rho^t = h + zeta_t(e) mod Z[i]
        if all(g_of(cmul(conj(c), cpow(RHO, t)))[0] == zeta(e, t) for t in range(-2, 3)):
            sel_zeta.add((a, b))
        dd = digits_at(c, [-1, 0, 1, 2])
        assert [dd[t] for t in (-1, 0, 1, 2)] == [digit7(e, t) for t in (-1, 0, 1, 2)]
    if all_values_in_2to5(c, G21):
        sel_vals.add((a, b))
print("type-7 points mod 325: text rule gives", len(sel_text), "points; values in {2..5}/7 on G(2,1):",
      sel_text == sel_vals, "; values at rho^t are zeta_t(3(a+bi)/7):", sel_zeta == sel_text)
print("  their digits at t=-1..2 are those of eps = 3(a+bi)/7 (aligned): True")
print("  points:", sorted(sel_text))
# compare with the SEVEN lines of cert_K2M1.txt (mod 325)
seven_pts = set()
for ln in open("run/cert_K2M1.txt"):
    if "SEVEN" in ln:
        x, y = ln.split("SEVEN")[1].split()
        X, Y = Fr(x), Fr(y)
        seven_pts.add((int((7 * X / 325) % 7), int((7 * Y / 325) % 7)))
print("  = the SEVEN points of cert_K2M1.txt:", seven_pts == sel_text)
# the same at modulus 65 (65 = 2 mod 7)
sel65 = sorted((a, b) for a, b in product(range(7), repeat=2) if ((2 * a) % 7, (2 * b) % 7) in A7)
print("type-7 points mod 65 by the same rule (2(a+bi) in A'_7):", sel65)
print("  matches the list after Proposition 8:",
      sel65 == sorted([(1, 2), (1, 5), (2, 1), (2, 6), (5, 1), (5, 6), (6, 2), (6, 5)]))
# lambda(a+bi) = 2a+3b and (65/7)(1+5i)
c = (Fr(65, 7), Fr(325, 7))
def red7(g):
    # residue of a rotation g (denominator prime to 7) in F_49
    den = 1
    for x in g:
        den = den * x.denominator // __import__('math').gcd(den, x.denominator)
    inv = pow(den, -1, 7)
    return ((int(g[0] * den) * inv) % 7, (int(g[1] * den) * inv) % 7)


G11 = [cmul(cpow(RHO, j), cpow(SIG, l)) for j in (-1, 0, 1) for l in (-1, 0, 1)]
print("(65/7)(1+5i): Re(conj(c) g) = lambda(g mod 7)/7 mod 1, lambda(a+bi)=2a+3b, for all g in G(1,1) (and i^a g):",
      all((cmul(conj(c), cmul(u, g))[0] - Fr((2 * red7(cmul(u, g))[0] + 3 * red7(cmul(u, g))[1]) % 7, 7)).denominator == 1
          for g in G11 for u in U4))

# the literal sentence "the character a -> Re(conj(eps) a) of Z_(7)[i] takes values in (1/7){2,..,5} + Z at every
# rational rotation"
eps = (Fr(2, 7), Fr(3, 7))
v = cmul(conj(eps), RHO)[0]
print("literal reading: eps=(2+3i)/7, a=rho: Re(conj(eps) rho) =", v, "-> in (1/7)Z?", (7 * v).denominator == 1)

# N = 5^k = 1 mod 7 iff 6 | k
print("5^k = 1 mod 7 exactly for 6|k (k<=36):", all((pow(5, k, 7) == 1) == (k % 6 == 0) for k in range(1, 37)))

# 2-adic points of the paragraph after Proposition 9
def least_margin(c, G):
    m = Fr(1)
    for g in G:
        w = cmul(conj(c), g)
        for x in (w[0], w[1]):
            f = x - fl(x)
            m = min(m, f, 1 - f)
    return m


for c in [(Fr(325, 4), Fr(325, 4)), (Fr(325, 4), Fr(325, 8)), (Fr(325, 8), Fr(325, 4))]:
    print("least margin on G(2,1) at", c, "=", least_margin(c, G21))

# index vectors of the 2-adic points vs the EXTRA vectors of cert_K2M1.txt
keys = None
vecs = {}
for ln in open("run/cert_K2M1.txt"):
    if ln.startswith("# functionals"):
        keys = [tuple(x.strip("()").split(",")) for x in ln.split(":", 1)[1].split()]
        keys = [(int(j), int(l), e) for (j, l, e) in keys]
    if ln.startswith("component"):
        head, lab = ln.split(" | ", 1)
        vecs[tuple(int(x) for x in head.split("indices")[1].split())] = lab.split()[0]


def ivec(c):
    out = []
    for (j, l, e) in keys:
        g = cmul(cpow(RHO, j), cpow(SIG, l))
        w = cmul(conj(c), g)
        out.append(fl(w[0]) if e == "Re" else fl(w[1]))
    return tuple(out)


for c in [(Fr(325, 4), Fr(325, 4)), (Fr(325, 4), Fr(325, 8)), (Fr(325, 8), Fr(325, 4)), (Fr(325, 8), Fr(325 * 3, 8))]:
    print("2-adic point", c, "least margin", least_margin(c, G21), "label of its index vector:", vecs.get(ivec(c)))
