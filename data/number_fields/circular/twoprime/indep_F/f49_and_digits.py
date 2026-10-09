"""Referee check of the finite facts (F1)-(F3) about F_49 = Z[i]/7, of the digit patterns of the 7-adic
characters, of the digits of the 13 type points modulo 325, and of the claim that any 3 consecutive digits
determine the family (c, q, 7) and the phase.  Exact arithmetic only (integers and Fractions).
Written from the definitions in the note; shares no code with the authors' programs."""
from fractions import Fraction as Fr
from itertools import product

P = 7
F49 = [(a, b) for a in range(P) for b in range(P)]


def mul(u, v):
    return ((u[0] * v[0] - u[1] * v[1]) % P, (u[0] * v[1] + u[1] * v[0]) % P)


def conj(u):
    return (u[0] % P, (-u[1]) % P)


def add(u, v):
    return ((u[0] + v[0]) % P, (u[1] + v[1]) % P)


mu8 = [z for z in F49 if mul(z, conj(z)) == (1, 0)]
print("mu8 =", mu8, "size", len(mu8))
assert len(mu8) == 8
Ap = [e for e in F49 if all(mul(conj(e), z)[0] in (2, 3, 4, 5) for z in mu8)]
print("A'_7 =", Ap, "size", len(Ap))
# (F1)
orbit = sorted(set(mul(Ap[0], z) for z in mu8))
print("(F1) single mu8-orbit:", orbit == sorted(Ap), "; 0 not in A'_7:", (0, 0) not in Ap)
claimed = sorted([(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)])
print("(F1) equals the list in the note:", sorted(Ap) == claimed)
# (F2): non-zero additive subgroups of F_49 = F_7-subspaces: 8 lines and the whole space
lines = set()
for d in F49:
    if d != (0, 0):
        lines.add(frozenset(((k * d[0]) % P, (k * d[1]) % P) for k in range(P)))
subgroups = list(lines) + [frozenset(F49)]
assert len(subgroups) == 9
bad = [(e, sorted(H)) for H in subgroups for e in F49 if all(add(e, h) in Ap for h in H)]
print("(F2) cosets of non-zero subgroups inside A'_7:", bad if bad else "none")
# (F3)
print("(F3) first coordinates of A'_7:", sorted(set(e[0] for e in Ap)))
print("conjugation-stable:", sorted(conj(e) for e in Ap) == sorted(Ap))
s = (0, 0)
for e in Ap:
    s = add(s, e)
print("sum of A'_7 =", s)

# rho, sigma mod 7
inv5 = pow(5, -1, 7)
inv13 = pow(13, -1, 7)
rho7 = ((3 * inv5) % 7, (4 * inv5) % 7)
sig7 = ((5 * inv13) % 7, (12 * inv13) % 7)
o = 1
z = rho7
while z != (1, 0):
    z = mul(z, rho7)
    o += 1
print("rho mod 7 =", rho7, " order", o, " rho^2 mod 7 =", mul(rho7, rho7), " sigma mod 7 =", sig7)
print("325 mod 7 =", 325 % 7, " 325^{-1} mod 7 =", pow(325, -1, 7), " 5^6 mod 7 =", pow(5, 6, 7),
      " order of 5 mod 7 =", [k for k in range(1, 7) if pow(5, k, 7) == 1][0])
pw = [(1, 0)]
for _ in range(7):
    pw.append(mul(pw[-1], rho7))
print("powers of rho mod 7 = mu8:", sorted(pw) == sorted(mu8))

# ---------------------------------------------------------------- digits
h = (Fr(1, 2), Fr(1, 2))
rho = (Fr(3, 5), Fr(4, 5))
rhoinv = (Fr(3, 5), Fr(-4, 5))


def cmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def csub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def cpow(u, j):
    r = (Fr(1), Fr(0))
    base = u if j >= 0 else (u[0], -u[1])  # |u| = 1
    for _ in range(abs(j)):
        r = cmul(r, base)
    return r


def floor(q):
    return q.numerator // q.denominator


def gvalue(w):
    """w = conj(c) rho^j.  Returns (g, n) with w = h + g + n, n in Z[i], g in [-1/2, 1/2)^2 (g = w - h - n)."""
    n = (floor(w[0] - h[0] + Fr(1, 2)), floor(w[1] - h[1] + Fr(1, 2)))
    g = (w[0] - h[0] - n[0], w[1] - h[1] - n[1])
    return g, n


def digits_of_point(c, ts):
    """digits nu_t = g_t - rho g_{t-1} of the complex number c (pair of Fractions) for t in ts."""
    cb = (c[0], -c[1])
    out = {}
    for t in ts:
        g1, _ = gvalue(cmul(cb, cpow(rho, t)))
        g0, _ = gvalue(cmul(cb, cpow(rho, t - 1)))
        out[t] = csub(g1, cmul(rho, g0))
    return out


def char_g(e, j):
    """g-value of the 7-adic character of class e (in F_49) at rho^j: (conj(e) rho~^j mod 7)/7 - h, reduced."""
    z = (1, 0)
    for _ in range(j % 8):
        z = mul(z, rho7)
    w = mul(conj(e), z)
    return (Fr(w[0], 7) - Fr(1, 2), Fr(w[1], 7) - Fr(1, 2))


def char_digit(e, j):
    return csub(char_g(e, j), cmul(rho, char_g(e, j - 1)))


gset = set()
pat7 = {}
for e in Ap:
    for j in range(8):
        gset.add(char_g(e, j))
    pat7[e] = [char_digit(e, j) for j in range(8)]
print("g-values of the 7-adic characters:", sorted(set(x for g in gset for x in g)))
L = [(Fr(2, 5), Fr(1, 5))]  # (2+i)/5
units = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))]
nonzero_ok = set(cmul(u, L[0]) for u in units)
S = pat7[Ap[0]]
print("S (x5) for e =", Ap[0], ":", [(int(5 * d[0]), int(5 * d[1])) for d in S])
ok_alt = all((S[j] == (0, 0)) != (S[(j + 1) % 8] == (0, 0)) for j in range(8))
ok_nz = all(d == (0, 0) or d in nonzero_ok for d in S)
ok_rot = all(S[(j + 2) % 8] == cmul((Fr(0), Fr(-1)), S[j]) for j in range(8))
print("zeros and non-zeros alternate:", ok_alt, "; non-zero digits are u(2+i)/5:", ok_nz, "; S_{j+2} = -i S_j:", ok_rot)
shifts = {}
for e in Ap:
    sh = [m for m in range(8) if all(pat7[e][j] == S[(j + m) % 8] for j in range(8))]
    shifts[e] = sh
print("each 7-pattern is a shift of S (shift amounts):", shifts)
print("the 8 shifts are distinct and exhaust Z/8:", sorted(x[0] for x in shifts.values()) == list(range(8)))

# q patterns: g_j = eta_j/6, eta_{j-1} = i eta_j  (Lemma values (iii)); digits nu_j = g_j - rho g_{j-1}
etas = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def q_g(eta0, j):
    # eta_j = (-i)^j eta_0
    e = (Fr(eta0[0]), Fr(eta0[1]))
    e = cmul(e, cpow((Fr(0), Fr(-1)), j % 4))
    return (e[0] / 6, e[1] / 6)


def q_digit(eta0, j):
    return csub(q_g(eta0, j), cmul(rho, q_g(eta0, j - 1)))


patq = {eta: [q_digit(eta, j) for j in range(8)] for eta in etas}
print("q digits all non-zero of the form u(2+i)/5:", all(d in nonzero_ok for p in patq.values() for d in p))
pats = {('c',): [(Fr(0), Fr(0))] * 8}
for eta in etas:
    pats[('q', eta)] = patq[eta]
for e in Ap:
    pats[('7', e)] = pat7[e]
assert len(pats) == 13
inj = True
for p in range(8):
    seen = {}
    for name, w in pats.items():
        tr = (w[p % 8], w[(p + 1) % 8], w[(p + 2) % 8])
        if tr in seen:
            inj = False
            print("  collision at position", p, name, seen[tr])
        seen[tr] = name
print("at every position, 3 consecutive digits determine the pattern among the 13 (c, 4 q, 8 seven):", inj)
# two digits do not suffice (sanity): check
inj2 = all(len(set((w[p % 8], w[(p + 1) % 8]) for w in pats.values())) == 13 for p in range(8))
print("(for comparison) 2 consecutive digits suffice:", inj2)

# ---------------------------------------------------------------- the 13 type points modulo 325
N0 = 325
types = [('c', (Fr(N0, 2), Fr(N0, 2)))]
for al in (1, 2):
    for be in (1, 2):
        types.append((('q', al, be), (Fr(N0 * al, 3), Fr(N0 * be, 3))))
for e in Ap:
    a, b = (5 * e[0]) % 7, (5 * e[1]) % 7
    types.append((('7', e), (Fr(N0 * a, 7), Fr(N0 * b, 7))))
# all rotations with denominator dividing 325
rots = [(Fr(x, N0), Fr(y, N0)) for x in range(-N0, N0 + 1) for y in range(-N0, N0 + 1) if x * x + y * y == N0 * N0]
print("rotations with denominator dividing 325:", len(rots))


def frac(q):
    return q - floor(q)


allpatterns_ok = True
for name, T in types:
    Tb = (T[0], -T[1])
    vals = set(frac(cmul(Tb, g)[0]) for g in rots)
    mn = min(min(v, 1 - v) for v in vals)
    dg = digits_of_point(T, [-1, 0, 1, 2])
    word = [dg[t] for t in (-1, 0, 1, 2)]
    # which bi-infinite pattern and phase has these 4 digits at positions -1..2 ?
    hits = [(nm, m) for nm, w in pats.items() for m in range(8) if all(w[(t + m) % 8] == dg[t] for t in (-1, 0, 1, 2))]
    if name[0] == '7':
        # expected: the character of class e itself (values at rho^t = character at rho~^t)
        exp_ok = all(dg[t] == char_digit(name[1], t) for t in (-1, 0, 1, 2))
        # its values at all 60 rotations are those of the character e (mod 1)
        chv = True
        for g in rots:
            # gamma mod 7
            gx = (g[0] * N0, g[1] * N0)
            gm = ((int(gx[0]) * pow(N0, -1, 7)) % 7, (int(gx[1]) * pow(N0, -1, 7)) % 7)
            w = mul(conj(name[1]), gm)
            if frac(cmul(Tb, g)[0]) != Fr(w[0], 7):
                chv = False
        extra = f" digits = those of character e: {exp_ok}; values = character e at all 60 rotations: {chv}"
        allpatterns_ok &= exp_ok and chv
    else:
        extra = ""
    print(f"type {name}: point/325 = ({T[0]/N0}, {T[1]/N0}); min margin over the 60 rotations = {mn};"
          f" values mod 1 = {sorted(vals)}; matching (pattern, shift): {hits[:3]}{extra}")
print("7-type points consistent with the characters:", allpatterns_ok)
