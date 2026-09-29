"""flat.py -- exact model of F = Q(zeta21) and the 126 unit directions U = D u conj(D).

Elements are integer 12-vectors v meaning (v_0 + v_1 z + ... + v_11 z^11)/7, z = zeta21,
Phi21(x) = x^12 - x^11 + x^9 - x^8 + x^6 - x^4 + x^3 - x + 1.
D = mu42 u omega*mu42 (Haugland), conj(D) = mu42 u conj(omega)*mu42, U = D u conj(D) (126 vectors).
Graph vertices are sums of directions, hence integer vectors (scaled by 7); all graph arithmetic is integral.
"""
import numpy as np
import cmath, math

N = 12
PHI = [1, -1, 0, 1, -1, 0, 1, 0, -1, 1, 0, -1, 1]      # coefficient of x^k at index k, degree 12, monic
OM7 = [4, -3, 1, 6, 1, -5, 1, 2, -7, 4, 1, -2]           # 7*omega in the power basis of zeta21


def red(r):
    """reduce an integer coefficient list (any length) mod Phi21 -> length-12 list"""
    r = list(r) + [0] * max(0, N - len(r))
    for d in range(len(r) - 1, N - 1, -1):
        c = r[d]
        if c:
            for k in range(N + 1):
                r[d - N + k] -= c * PHI[k]
            assert r[d] == 0
    return r[:N]


def mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    return red(r)


def zpow(k):
    e = [0] * (k % 21 + 1)
    e[k % 21] = 1
    return red(e)


# conjugation z -> z^20 = z^-1 as an integer matrix acting on coefficient vectors
CONJ = np.array([zpow(-k) for k in range(N)], dtype=np.int64).T   # column k = image of z^k


def conj(v):
    return tuple(int(x) for x in CONJ @ np.array(v, dtype=np.int64))


Z21 = cmath.exp(2j * math.pi / 21)
ZP = np.array([Z21 ** k for k in range(N)])


def cval(v):
    """complex value of the element v/7 under z -> exp(2 pi i/21)"""
    return complex(np.dot(np.asarray(v, dtype=float), ZP) / 7.0)


def directions():
    """return (D, CD, U) as lists of integer 12-tuples (scaled by 7).
    D[j] = u_j of Haugland: u_{2m} = zeta42^m, u_{2m+1} = zeta42^m omega; zeta42 = -z^11.
    CD[j] = conj(u_j).  U = D followed by the 42 vectors conj(omega)*mu42 (dedup, order kept)."""
    D = []
    for j in range(84):
        m = j // 2
        base = [(-1) ** m * c for c in zpow(11 * m)]
        if j % 2 == 0:
            D.append(tuple(7 * c for c in base))
        else:
            D.append(tuple(mul(base, OM7)))
    CD = [conj(u) for u in D]
    U, seen = [], set()
    for u in D + CD:
        if u not in seen:
            seen.add(u)
            U.append(u)
    return D, CD, U


# ---------------- hashing of integer vectors (exact lookup via verified random linear keys) -------------
_rng = np.random.default_rng(12345)
HW = _rng.integers(1, 2 ** 62, size=N, dtype=np.int64)


def keys(A):
    A = np.asarray(A, dtype=np.int64)
    with np.errstate(over='ignore'):
        return (A * HW).sum(axis=1)          # wraps mod 2^64; exactness is checked on every hit


class PointSet:
    """a set of integer 12-vectors with exact membership queries"""

    def __init__(self, P):
        self.P = np.asarray(P, dtype=np.int64).reshape(-1, N)
        k = keys(self.P)
        order = np.argsort(k, kind='stable')
        self.sk = k[order]
        self.order = order
        if len(k) > 1 and np.any(self.sk[1:] == self.sk[:-1]):
            # either duplicate points or a hash collision; both are errors for our use
            raise ValueError("duplicate keys among points")

    def index(self, Q):
        """for each row of Q: index into P or -1 (exact)"""
        Q = np.asarray(Q, dtype=np.int64).reshape(-1, N)
        k = keys(Q)
        pos = np.searchsorted(self.sk, k)
        pos[pos >= len(self.sk)] = 0
        hit = self.sk[pos] == k if len(self.sk) else np.zeros(len(k), bool)
        idx = np.where(hit, self.order[pos] if len(self.sk) else -1, -1)
        ok = idx >= 0
        if ok.any():
            same = np.all(self.P[idx[ok]] == Q[ok], axis=1)
            tmp = idx[ok]
            tmp[~same] = -1
            idx[ok] = tmp
        return idx


def unique_rows(A):
    A = np.asarray(A, dtype=np.int64).reshape(-1, N)
    return np.unique(A, axis=0)


def build_edges(P, U):
    """all pairs (a<b) with P[b]-P[a] in U; returns int array (m,2) and array of direction indices"""
    P = np.asarray(P, dtype=np.int64)
    S = PointSet(P)
    E, J = [], []
    Ua = np.asarray(U, dtype=np.int64)
    for j, u in enumerate(Ua):
        idx = S.index(P + u)
        a = np.nonzero(idx >= 0)[0]
        b = idx[a]
        keep = a < b
        E.append(np.stack([a[keep], b[keep]], axis=1))
        J.append(np.full(keep.sum(), j))
    return np.concatenate(E), np.concatenate(J)


def sums_ball(dirs, k, extra=()):
    """sums of at most k vectors of dirs (as a unique int array), including 0"""
    Da = np.asarray(dirs, dtype=np.int64)
    cur = np.zeros((1, N), dtype=np.int64)
    allp = [cur]
    for _ in range(k):
        nxt = (cur[:, None, :] + Da[None, :, :]).reshape(-1, N)
        nxt = unique_rows(nxt)
        cur = nxt
        allp.append(cur)
    return unique_rows(np.concatenate(allp))


def conj_rows(A):
    A = np.asarray(A, dtype=np.int64)
    return (CONJ @ A.T).T.copy()


def find_triangle(n, E, prefer=None):
    """find a triangle (a,b,c) in the graph; prefer containing vertex `prefer`"""
    adj = [set() for _ in range(n)]
    for a, b in E:
        adj[a].add(b); adj[b].add(a)
    order = [prefer] if prefer is not None else []
    order += list(range(n))
    for a in order:
        for b in sorted(adj[a]):
            common = adj[a] & adj[b]
            if common:
                return (a, b, min(common))
    return None


def cnf_clauses(n, E, tri, k=4, selectors=False):
    """4-colouring CNF: var(v,c) = 1 + k*v + c. At-least-one per vertex, one conflict clause per edge
    and colour, the colours of triangle tri fixed to 0,1,2. If selectors: selector of v = 1 + k*n + v,
    ALO clause becomes (-s_v | ...)."""
    cls = []
    for v in range(n):
        c = [1 + k * v + i for i in range(k)]
        if selectors:
            c = [-(1 + k * n + v)] + c
        cls.append(c)
    for a, b in E:
        for i in range(k):
            cls.append([-(1 + k * a + i), -(1 + k * b + i)])
    if tri is not None:
        for i, v in enumerate(tri):
            cls.append([1 + k * v + i])
    return cls
