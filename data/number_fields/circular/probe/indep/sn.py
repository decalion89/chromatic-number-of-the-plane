"""S_N^(r) mod N, N = 5^k, computed exactly by lifting (level k-1 pieces, 25 translates, cut at j = +-k), plus
independent pointwise tests (G_N enumerated by brute force over |alpha|^2 = N^2) and the claimed structure of
Proposition D1 (P_k^(s) around c*_N and Y_k(q) around q in Q_N), built from scratch as H-polygons."""
from fractions import Fraction as Fr
from math import isqrt
from hgeom import Poly, funcs, ev, fl, gmul, gconj, rho_pow


def square(r):
    lo, hi = r, 1 - r
    return Poly([(Fr(1), Fr(0), hi), (Fr(-1), Fr(0), -lo), (Fr(0), Fr(1), hi), (Fr(0), Fr(-1), -lo)])


def lift_levels(K, r, keep_all=True):
    lo, hi = r, 1 - r
    levels = [[square(r)]]
    for k in range(1, K + 1):
        Np = 5 ** (k - 1)
        new = []
        for P in levels[-1]:
            for a in range(5):
                for b in range(5):
                    pieces = [P.translate((Fr(Np * a), Fr(Np * b)))]
                    for j in (k, -k):
                        for f in funcs(j):
                            nxt = []
                            for R in pieces:
                                nxt.extend(R.strip_split(f, lo, hi))
                            pieces = nxt
                            if not pieces:
                                break
                        if not pieces:
                            break
                    new.extend(pieces)
        levels.append(new)
        if not keep_all:
            levels[-2] = None
    return levels


def direct_level(k, r):
    """independent second enumeration: all N^2 unit squares of S_1 in [0,N)^2 cut by every |j| <= k at once"""
    N = 5 ** k
    lo, hi = r, 1 - r
    sq = square(r)
    out = []
    for a in range(N):
        for b in range(N):
            pieces = [sq.translate((Fr(a), Fr(b)))]
            for j in range(1, k + 1):
                for jj in (j, -j):
                    for f in funcs(jj):
                        nxt = []
                        for R in pieces:
                            nxt.extend(R.strip_split(f, lo, hi))
                        pieces = nxt
                        if not pieces:
                            break
                    if not pieces:
                        break
                if not pieces:
                    break
            out.extend(pieces)
    return out


_GN = {}


def GN_brute(N):
    if N in _GN:
        return _GN[N]
    out = []
    for a in range(-N, N + 1):
        b2 = N * N - a * a
        b = isqrt(b2)
        if b * b == b2:
            for bb in sorted({b, -b}):
                out.append((Fr(a, N), Fr(bb, N)))
    _GN[N] = out
    return out


def GN_rho(k):
    out = []
    for j in range(-k, k + 1):
        z = rho_pow(j)
        for e in ((1, 0), (0, 1), (-1, 0), (0, -1)):
            out.append(gmul((Fr(e[0]), Fr(e[1])), z))
    return out


def dist_int(t):
    f = t - fl(t)
    return min(f, 1 - f)


def kappa(c, G):
    """min over gamma of || Re(conj(c) gamma) ||"""
    return min(dist_int(c[0] * g[0] + c[1] * g[1]) for g in G)


def in_S(c, G, r):
    return kappa(c, G) >= r


def cstar(N):
    return (Fr(N, 2), Fr(N, 2))


def QN(N):
    return [(Fr(N * a, 3), Fr(N * b, 3)) for a in (1, 2) for b in (1, 2)]


def Pk(k, s):
    H = []
    for j in range(-k, k + 1):
        for f in funcs(j):
            H.append((f[0], f[1], s))
            H.append((-f[0], -f[1], s))
    return Poly(H)


def Yk(q, k, r):
    """y with conj(q+y) rho^j in the same lift square [n+r, n+1-r]^2 as conj(q) rho^j, |j| <= k"""
    H = []
    for j in range(-k, k + 1):
        for f in funcs(j):
            val = ev(f, q)
            n = fl(val)
            fr = val - n
            assert fr in (Fr(1, 3), Fr(2, 3)), (q, j, val)
            # n + r <= val + f(y) <= n + 1 - r
            H.append((f[0], f[1], n + 1 - r - val))
            H.append((-f[0], -f[1], val - n - r))
    return Poly(H)


def nearest_shift(p, T, N):
    """lattice vector N*m with p - T - N*m smallest (componentwise rounding)"""
    m0 = round((p[0] - T[0]) / N)
    m1 = round((p[1] - T[1]) / N)
    return (T[0] + N * m0, T[1] + N * m1)


def classify(P, k, r, refs=None):
    """return ('c', None) / ('q', q) if P equals c*_N + P_k^(s) / q + Y_k(q) modulo N, else None"""
    N = 5 ** k
    s = Fr(1, 2) - r
    if refs is None:
        refs = make_refs(k, r)
    cen = P.centroid_v()
    for lab, T, ref in refs:
        Tm = nearest_shift(cen, T, N)
        sh = sorted((v[0] - Tm[0], v[1] - Tm[1]) for v in P.V)
        if sh == ref:
            return (lab, T)
    return None


def make_refs(k, r):
    N = 5 ** k
    s = Fr(1, 2) - r
    refs = [("c", cstar(N), sorted(Pk(k, s).V))]
    for q in QN(N):
        refs.append(("q", q, sorted(Yk(q, k, r).V)))
    return refs
