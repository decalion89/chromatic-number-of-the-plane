"""Square roots in a multiquadratic field, by the obvious recursion.

F = Q(sqrt(g_1), ..., sqrt(g_k)) with the basis indexed by subsets, so an
element splits as y = a + b.sqrt(g_k) with a, b in the field one generator
shorter.  If sqrt(y) = u + v.sqrt(g_k) then a = u^2 + g_k v^2 and b = 2uv, so
u^2 and g_k v^2 are the roots of

    T^2 - a T + g_k b^2 / 4,

whose discriminant a^2 - g_k b^2 must itself be a square one level down.  The
base case is Q.  The b = 0 branch needs care: y may be a square, or g_k times
one, and both have to be tried.

Everything is exact; nothing is reconstructed from floats.
"""
from fractions import Fraction as Fr
from math import isqrt


def madd(x, y):
    return tuple(a + b for a, b in zip(x, y))


def msub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def mscal(q, x):
    q = Fr(q)
    return tuple(q * a for a in x)


def mmul(x, y, gens):
    d = len(x)
    out = [Fr(0)] * d
    for i in range(d):
        if not x[i]:
            continue
        for j in range(d):
            if not y[j]:
                continue
            common, fac = i & j, 1
            for b in range(len(gens)):
                if common >> b & 1:
                    fac *= gens[b]
            out[i ^ j] += x[i] * y[j] * fac
    return tuple(out)


def minv(x, gens):
    d = len(x)
    cols = []
    for i in range(d):
        e = [Fr(0)] * d
        e[i] = Fr(1)
        cols.append(mmul(x, tuple(e), gens))
    A = [[cols[j][i] for j in range(d)] + [Fr(1 if i == 0 else 0)]
         for i in range(d)]
    for c in range(d):
        p = next(t for t in range(c, d) if A[t][c])
        A[c], A[p] = A[p], A[c]
        sc = Fr(1) / A[c][c]
        A[c] = [v * sc for v in A[c]]
        for t in range(d):
            if t != c and A[t][c]:
                f = A[t][c]
                A[t] = [u - f * v for u, v in zip(A[t], A[c])]
    return tuple(A[i][d] for i in range(d))


def _ratsqrt(q):
    if q < 0:
        return None
    n, d = q.numerator, q.denominator
    rn, rd = isqrt(n), isqrt(d)
    return Fr(rn, rd) if rn * rn == n and rd * rd == d else None


def msqrt(c, gens):
    """A square root of c in Q(sqrt g for g in gens), or None."""
    k = len(gens)
    if k == 0:
        r = _ratsqrt(c[0])
        return (r,) if r is not None else None
    half = 1 << (k - 1)
    a, b, g, sub = c[:half], c[half:], gens[-1], gens[:-1]
    zero = (Fr(0),) * half
    if all(x == 0 for x in b):
        r = msqrt(a, sub)
        if r is not None:
            return r + zero
        r = msqrt(mscal(Fr(1, g), a), sub)     # y = g . square
        if r is not None:
            return zero + r
        return None
    delta = msub(mmul(a, a, sub), mscal(g, mmul(b, b, sub)))
    s = msqrt(delta, sub)
    if s is None:
        return None
    for t in (mscal(Fr(1, 2), madd(a, s)), mscal(Fr(1, 2), msub(a, s))):
        u = msqrt(t, sub)
        if u is None or all(x == 0 for x in u):
            continue
        v = mmul(b, minv(mscal(2, u), sub), sub)
        cand = u + v
        if mmul(cand, cand, gens) == c:
            return cand
    return None
