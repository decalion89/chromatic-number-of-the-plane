"""Digit dynamics of the probe (exact arithmetic).

A character chi of A = Z[1/5][i] is the same as a sequence of values V_j = conj(c) rho^j mod Z[i] (j in Z) with
V_{j+1} = rho V_j mod Lambda_+ (Lambda_+ = ((2+i)/5) Z[i]); the condition kappa(chi) >= theta says
V_j in h + B_s + Z[i] for all j, where s = 1/2 - theta, B_s = [-s, s]^2, h = (1+i)/2.  Writing V_j = h + g_j with
g_j in B_s (unique for s < 1/2), we get

        g_{j+1} = rho g_j + nu_{j+1},     nu_{j+1} in Lambda_+,  g_j in B_s.

For s <= 3/14 the digit nu has sup-norm <= s(1 + sqrt2) < 0.52, so nu is one of five 'digits'
        '0' = 0,  'a' = (2+i)/5,  'b' = i(2+i)/5,  'c' = -(2+i)/5,  'd' = -i(2+i)/5.
The probe S_N^theta (N = 5^k) is the set of finite sequences (g_j)_{|j|<=k}; its connected components correspond
one-to-one to digit words (nu_{-k+1}, ..., nu_k) whose polygon
        P(word) = { g_0 : g_j in B_s for |j| <= k }
is non-empty (each g_j is an affine function rho^j g_0 + c_j of g_0).  The main components are: type c (all digits
'0'; the point N h + x) and type q (the pattern nu_{j+1} = -i nu_j, i.e. a -> d -> c -> b -> a; the four points of
Q_N).  Everything here is exact (fractions.Fraction).
"""
from fractions import Fraction as Fr

ZERO = Fr(0)
DIG = {'0': (Fr(0), Fr(0)), 'a': (Fr(2, 5), Fr(1, 5)), 'b': (Fr(-1, 5), Fr(2, 5)),
       'c': (Fr(-2, 5), Fr(-1, 5)), 'd': (Fr(1, 5), Fr(-2, 5))}
LETTERS = '0abcd'
QNEXT = {'a': 'd', 'd': 'c', 'c': 'b', 'b': 'a'}      # nu_{j+1} = -i nu_j
QPREV = {v: k for k, v in QNEXT.items()}


def rho(p):
    x, y = p
    return ((3 * x - 4 * y) / 5, (4 * x + 3 * y) / 5)


def rhoinv(p):
    x, y = p
    return ((3 * x + 4 * y) / 5, (-4 * x + 3 * y) / 5)


def add(p, q):
    return (p[0] + q[0], p[1] + q[1])


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1])


def cmul(p, q):
    return (p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0])


def rho_pow(j):
    z = (Fr(1), Fr(0))
    for _ in range(abs(j)):
        z = rho(z) if j > 0 else rhoinv(z)
    return z


# ------------------------------------------------------------------ exact convex clipping
def clip(poly, a, b, c):
    """part of the convex polygon (vertex list) with a x + b y <= c (exact; keeps degenerate pieces)."""
    out = []
    n = len(poly)
    if n == 0:
        return []
    vals = [c - (a * p[0] + b * p[1]) for p in poly]
    if n == 1:
        return list(poly) if vals[0] >= 0 else []
    for i in range(n):
        P, Q = poly[i], poly[(i + 1) % n]
        fp, fq = vals[i], vals[(i + 1) % n]
        if fp >= 0:
            out.append(P)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    res = []
    for p in out:
        if not res or res[-1] != p:
            res.append(p)
    while len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res


def clip_square(poly, A, B, cvec, s):
    """keep g_0 in poly with g = (A x - B y + c1, B x + A y + c2) in [-s, s]^2"""
    c1, c2 = cvec
    for (a, b, c) in ((A, -B, s - c1), (-A, B, s + c1), (B, A, s - c2), (-B, -A, s + c2)):
        poly = clip(poly, a, b, c)
        if not poly:
            return []
    return poly


def square(s):
    return [(-s, -s), (s, -s), (s, s), (-s, s)]


def word_type(word):
    """'C' (all zero), 'Q' (q pattern throughout, all non-zero) or 'X'."""
    if all(ch == '0' for ch in word):
        return 'C'
    if all(ch != '0' for ch in word) and all(QNEXT[word[i]] == word[i + 1] for i in range(len(word) - 1)):
        return 'Q'
    return 'X'


def runs(word):
    """decompose a word into maximal runs: ('0', length) for zero runs and ('q', length) for q-pattern runs
    (a maximal stretch of non-zero digits with nu_{j+1} = -i nu_j); a break of the q-pattern between two non-zero
    digits starts a new q-run (marked 'q!')."""
    out = []
    i = 0
    n = len(word)
    while i < n:
        if word[i] == '0':
            j = i
            while j < n and word[j] == '0':
                j += 1
            out.append(('0', j - i))
            i = j
        else:
            j = i + 1
            while j < n and word[j] != '0' and QNEXT[word[j - 1]] == word[j]:
                j += 1
            tag = 'q'
            if out and out[-1][0] in ('q', 'q!'):
                tag = 'q!'
            out.append((tag, j - i))
            i = j
    return out


# ------------------------------------------------------------------ the probe as words (two-sided lift)
class Comp:
    __slots__ = ('word', 'poly', 'cR', 'cL', 'k')

    def __init__(self, word, poly, cR, cL, k):
        self.word, self.poly, self.cR, self.cL, self.k = word, poly, cR, cL, k


def level0(s):
    return [Comp('', square(s), (Fr(0), Fr(0)), (Fr(0), Fr(0)), 0)]


def lift(comps, s):
    """from level k-1 (window [-(k-1), k-1], digits nu_{-k+2..k-1}) to level k: add nu_k on the right and
    nu_{-k+1} on the left."""
    out = []
    if not comps:
        return out
    k = comps[0].k + 1
    AR, BR = rho_pow(k)
    AL, BL = rho_pow(-k)
    for C in comps:
        for r in LETTERS:
            cR = add(rho(C.cR), DIG[r])                 # g_k = rho g_{k-1} + nu_k
            pR = clip_square(C.poly, AR, BR, cR, s)
            if not pR:
                continue
            for l in LETTERS:
                cL = rhoinv(sub(C.cL, DIG[l]))           # g_{-k} = rho^{-1}(g_{-k+1} - nu_{-k+1})
                p = clip_square(pR, AL, BL, cL, s)
                if p:
                    out.append(Comp(l + C.word + r, p, cR, cL, k))
    return out


def forward_lift(comps, s):
    """one-sided: add a digit on the right only (window [k0, k])."""
    out = []
    if not comps:
        return out
    k = comps[0].k + 1
    AR, BR = rho_pow(k)
    for C in comps:
        for r in LETTERS:
            cR = add(rho(C.cR), DIG[r])
            pR = clip_square(C.poly, AR, BR, cR, s)
            if pR:
                out.append(Comp(C.word + r, pR, cR, C.cL, k))
    return out


def word_polygon(word, s, j0=0):
    """polygon of g_{j0} for a word nu_{j0+1} ... nu_{j0+n} (window [j0, j0+n]), in g_{j0} coordinates"""
    poly = square(s)
    c = (Fr(0), Fr(0))
    for t, ch in enumerate(word, start=1):
        c = add(rho(c), DIG[ch])
        A, B = rho_pow(t)
        poly = clip_square(poly, A, B, c, s)
        if not poly:
            return []
    return poly


def orbit(g0, word):
    """the points g_0, g_1, ..., g_n for a word nu_1..nu_n"""
    pts = [g0]
    for ch in word:
        pts.append(add(rho(pts[-1]), DIG[ch]))
    return pts


def supnorm(p):
    return max(abs(p[0]), abs(p[1]))
