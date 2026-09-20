"""Build the chain over a field chosen to block, and measure it.

The tower is K = L(W)(Y) with

    L = Q[T]/(T^3 - 9T^2 + 14T + 8),   m = T,
    W^2 = V(m) = -39m^4 + 540m^3 - 666m^2 - 324m + 297,
    Y^2 = -3,

so F = L(W) is totally real of degree 6 -- V is positive at all three
embeddings -- and K = F(Y) is its CM quadratic extension, degree 12, with
complex conjugation Y -> -Y.  The cubic was chosen so that every prime of F
above 5 has residue degree at least 3, which is exactly the condition the
residue-degree theorem requires of a graph that blocks.

Inside it the chain is explicit:

    r  = (m^2 - 6m - 3)/(m^2 + 3)        |r| <= 2 automatically
    y1 = 3 + m(r - 1)                    y1^2 = 3(4 - r^2)
    y2 = W/(m^2 + 3)                     y2^2 = 3(8 - 12r - 9r^2)
    a  = (r + y1 Y/3)/2                  |a| = 1
    b  = (-8/3 - r + y2 Y/9) / (2(1 + abar))
    zeta_6 = (1 + Y)/2

and |1 + a + b|^2 = 1/3, so the three rhombi on 1, a, b close.
"""
import sys, time, cmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

CUB = (8, 14, -9)          # T^3 - 9T^2 + 14T + 8  (ascending constant first)


def l_zero():
    return (Fr(0),) * 3


def l_rat(q):
    return (Fr(q), Fr(0), Fr(0))


L_T = (Fr(0), Fr(1), Fr(0))


def l_add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def l_sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def l_scal(c, x):
    return tuple(Fr(c) * a for a in x)


def l_mul(x, y):
    r = [Fr(0)] * 5
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                r[i + j] += a * b
    for k in (4, 3):
        c = r[k]
        if c:
            r[k] = Fr(0)
            r[k - 3] -= c * CUB[0]
            r[k - 2] -= c * CUB[1]
            r[k - 1] -= c * CUB[2]
    return tuple(r[:3])


def l_inv(x):
    M = [[l_mul(x, (Fr(1) if k == j else Fr(0),
                    Fr(1) if k == j and False else Fr(0),
                    Fr(0)))[0] for k in range(3)] for j in range(3)]
    cols = []
    for j in range(3):
        e = tuple(Fr(1) if k == j else Fr(0) for k in range(3))
        cols.append(l_mul(x, e))
    A = [[cols[j][i] for j in range(3)] + [Fr(1) if i == 0 else Fr(0)]
         for i in range(3)]
    for c in range(3):
        p = next(r for r in range(c, 3) if A[r][c])
        A[c], A[p] = A[p], A[c]
        s = Fr(1) / A[c][c]
        A[c] = [v * s for v in A[c]]
        for r in range(3):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [u - f * v for u, v in zip(A[r], A[c])]
    return tuple(A[i][3] for i in range(3))


MM = l_mul(L_T, L_T)
V = l_add(l_add(l_scal(-39, l_mul(MM, MM)), l_scal(540, l_mul(MM, L_T))),
          l_add(l_scal(-666, MM), l_add(l_scal(-324, L_T), l_rat(297))))


def m_of(x):
    return (x, l_zero())


def m_add(x, y):
    return (l_add(x[0], y[0]), l_add(x[1], y[1]))


def m_sub(x, y):
    return (l_sub(x[0], y[0]), l_sub(x[1], y[1]))


def m_mul(x, y):
    return (l_add(l_mul(x[0], y[0]), l_mul(l_mul(x[1], y[1]), V)),
            l_add(l_mul(x[0], y[1]), l_mul(x[1], y[0])))


def m_inv(x):
    d = l_sub(l_mul(x[0], x[0]), l_mul(l_mul(x[1], x[1]), V))
    di = l_inv(d)
    return (l_mul(x[0], di), l_mul(l_scal(-1, x[1]), di))


M_ONE = m_of(l_rat(1))
M_W = (l_zero(), l_rat(1))


def k_of(x):
    return (x, (l_zero(), l_zero()))


def k_add(x, y):
    return (m_add(x[0], y[0]), m_add(x[1], y[1]))


def k_sub(x, y):
    return (m_sub(x[0], y[0]), m_sub(x[1], y[1]))


def k_mul(x, y):
    ac = m_mul(x[0], y[0])
    bd = m_mul(x[1], y[1])
    return (m_sub(ac, (l_scal(3, bd[0]), l_scal(3, bd[1]))),
            m_add(m_mul(x[0], y[1]), m_mul(x[1], y[0])))


def k_conj(x):
    return (x[0], (l_scal(-1, x[1][0]), l_scal(-1, x[1][1])))


def k_inv(x):
    d = m_add(m_mul(x[0], x[0]), (l_scal(3, m_mul(x[1], x[1])[0]),
                                  l_scal(3, m_mul(x[1], x[1])[1])))
    di = m_inv(d)
    return (m_mul(x[0], di), m_mul((l_scal(-1, x[1][0]),
                                    l_scal(-1, x[1][1])), di))


def k_norm2(x):
    return k_mul(x, k_conj(x))[0]


def k_rat(q):
    return k_of(m_of(l_rat(q)))


K_ONE = k_rat(1)
K_ZERO = k_rat(0)
K_Y = (m_of(l_zero()), M_ONE)

mm = k_of(m_of(L_T))
mm2 = k_mul(mm, mm)
den = k_add(mm2, k_rat(3))
r = k_mul(k_sub(k_sub(mm2, k_mul(k_rat(6), mm)), k_rat(3)), k_inv(den))
y1 = k_add(k_rat(3), k_mul(mm, k_sub(r, K_ONE)))
assert k_mul(y1, y1) == k_mul(k_rat(3), k_sub(k_rat(4), k_mul(r, r)))
y2 = k_mul(k_of(M_W), k_inv(den))
lhs = k_mul(y2, y2)
rhs = k_mul(k_rat(3), k_sub(k_sub(k_rat(8), k_mul(k_rat(12), r)),
                            k_mul(k_rat(9), k_mul(r, r))))
assert lhs == rhs
HALF = k_rat(Fr(1, 2))
a = k_mul(HALF, k_add(r, k_mul(y1, k_mul(k_rat(Fr(1, 3)), K_Y))))
assert k_norm2(a) == M_ONE, "a must be a unit step"
Tv = k_sub(k_rat(Fr(-8, 3)), r)
b = k_mul(k_add(Tv, k_mul(y2, k_mul(k_rat(Fr(1, 9)), K_Y))),
          k_inv(k_mul(k_rat(2), k_add(K_ONE, k_conj(a)))))
print(f"|a|^2 = 1 ? {k_norm2(a) == M_ONE};  |b|^2 = 1 ? "
      f"{k_norm2(b) == M_ONE}", flush=True)
s3 = k_add(k_add(K_ONE, a), b)
print(f"|1 + a + b|^2 = 1/3 ? "
      f"{k_norm2(s3) == m_of(l_rat(Fr(1, 3)))}", flush=True)
Z6 = k_mul(HALF, k_add(K_ONE, K_Y))
assert k_norm2(Z6) == M_ONE
ONE_PLUS = k_add(K_ONE, Z6)
assert k_norm2(ONE_PLUS) == m_of(l_rat(3))

pts, B = [K_ZERO], K_ZERO
for w in (K_ONE, a, b):
    for off in (w, k_mul(Z6, w), k_mul(ONE_PLUS, w)):
        pts.append(k_add(B, off))
    B = k_add(B, k_mul(ONE_PLUS, w))
uniq, seen = [], set()
for q in pts:
    if q not in seen:
        seen.add(q)
        uniq.append(q)
E = [(i, j) for i in range(len(uniq)) for j in range(i + 1, len(uniq))
     if k_norm2(k_sub(uniq[j], uniq[i])) == M_ONE]
chi = None
for kk in (3, 4, 5):
    cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(uniq))]
    for x, y in E:
        for c in range(kk):
            cls.append([-(1 + x * kk + c), -(1 + y * kk + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if sv.solve():
            chi = kk
            break


def flat(x):
    out = []
    for mpart in x:
        for lpart in mpart:
            out.extend(lpart)
    return tuple(out)


vecs = []
for x, y in E:
    vecs.append(flat(k_sub(uniq[y], uniq[x])))
    vecs.append(flat(k_sub(uniq[x], uniq[y])))
dn = 1
for v in vecs:
    for q in v:
        dn = dn * q.denominator // gcd(dn, q.denominator)
iv = set()
for v in vecs:
    w = tuple(int(q * dn) for q in v)
    g = 0
    for t in w:
        g = gcd(g, abs(t))
    iv.add(tuple(t // g for t in w) if g > 1 else w)
rows = [[Fr(x) for x in v] for v in iv]
rank = 0
for c in range(12):
    pr = next((t for t in range(rank, len(rows)) if rows[t][c]), None)
    if pr is None:
        continue
    rows[rank], rows[pr] = rows[pr], rows[rank]
    f = rows[rank][c]
    rows[rank] = [x / f for x in rows[rank]]
    for t in range(len(rows)):
        if t != rank and rows[t][c]:
            k2 = rows[t][c]
            rows[t] = [x - k2 * y for x, y in zip(rows[t], rows[rank])]
    rank += 1
phi, _ = has_homomorphism(sorted(iv), 5)
print(f"chain over the degree-12 tower: {len(uniq)} points, {len(E)} edges, "
      f"chi = {chi}, {len(iv)} directions, module rank {rank}", flush=True)
print("  " + ("*** BLOCKS ***" if phi is None else "coset colouring exists"),
      flush=True)
