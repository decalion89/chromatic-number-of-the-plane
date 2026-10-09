"""limit.py: the N -> infinity limit of the configuration G_N * Delta' (Delta' = {1} u Delta), N = 5^k.
A unit vector of K^2 = Q(i) + sqrt(d) Q(i) is v = X + sqrt(d) Y.  A character is (X, Y) -> Re(conj(c1) X) + Re(conj(c2) Y);
on the slice G_N * delta it is gamma -> Re(conj(w_delta) gamma) with w_delta = c1 conj(delta_X) + c2 conj(delta_Y).
By the structure lemma S_N = (N(1+i)/2 + P) u Q mod N with P of bounded size, so for large N the configuration is
feasible iff some C in Q(i)^2 has w_delta(C) in E + Z[i] for every delta in Delta', E = {(1+i)/2} u {(a+bi)/3: a,b in {1,2}}.
With M the real 2n x 4 matrix of C -> (w_delta(C)), the solutions mod Z^{2n} of 6 w in Z^{2n} form L/6L where
L = im(M) cap Z^{2n}; we enumerate L mod 6 (6^4 elements) and test the pattern."""
import sys, itertools
from fractions import Fraction as F
from math import lcm
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ


def int_kernel(A, m):
    """integer basis of {x in Z^m : A x = 0} for an integer matrix A (list of rows), via LLL on [I | W A^T]"""
    W = 10 ** 8
    rows = [[1 if j == i else 0 for j in range(m)] + [W * A[r][i] for r in range(len(A))] for i in range(m)]
    M = DomainMatrix([[ZZ(x) for x in r] for r in rows], (m, m + len(A)), ZZ).lll()
    R = [list(map(int, M.to_Matrix().row(i))) for i in range(m)]
    ker = [r[:m] for r in R if all(x == 0 for x in r[m:])]
    for r in ker:
        assert all(sum(A[q][i] * r[i] for i in range(m)) == 0 for q in range(len(A)))
    return ker


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def conj(a):
    return (a[0], -a[1])


def matrix(deltas):
    """rows: Re and Im of w_delta as linear forms in C = (c1r, c1i, c2r, c2i)"""
    M = []
    for dX, dY in deltas:
        a, b = conj(dX), conj(dY)
        # C1 * a: re = c1r a0 - c1i a1, im = c1r a1 + c1i a0
        M.append([a[0], -a[1], b[0], -b[1]])
        M.append([a[1], a[0], b[1], b[0]])
    return M


def solutions(deltas, want_all=False):
    M = matrix(deltas); n2 = len(M)
    D = 1
    for r in M:
        for x in r: D = lcm(D, F(x).denominator)
    A = [[int(F(x) * D) for x in r] for r in M]
    At = [[A[r][c] for r in range(n2)] for c in range(4)]          # 4 x n2
    K = int_kernel(At, n2)                                         # left kernel of M
    L = int_kernel(K, n2) if K else [[1 if j == i else 0 for j in range(n2)] for i in range(n2)]
    assert len(L) == 4, len(L)
    sols = []
    for y in itertools.product(range(6), repeat=4):
        v = [sum(y[j] * L[j][i] for j in range(4)) % 6 for i in range(n2)]
        ok = True
        pat = []
        for q in range(0, n2, 2):
            p = (v[q], v[q + 1])
            if p == (3, 3): pat.append("c")
            elif p[0] in (2, 4) and p[1] in (2, 4): pat.append("q")
            else: ok = False; break
        if ok:
            sols.append("".join(pat))
            if not want_all: return sols
    return sols


def ntype(d, n):
    """u_n = ((d - n^2) + 2 i n sqrt d)/(d + n^2) as (delta_X, delta_Y)"""
    t = F(d + n * n)
    return ((F(d - n * n) / t, F(0)), (F(0), F(2 * n) / t))


def utype(d, n, s=1):
    """u_n = (n + i sqrt d)^2/(n^2 + d) = ((n^2 - d) + 2 i n sqrt d)/(n^2 + d) (s = -1: its mirror image)"""
    t = F(n * n + d)
    return ((F(n * n - d) / t, F(0)), (F(0), F(2 * n * s) / t))


def tconj(u):
    (x, y) = u
    return (x, (-y[0], -y[1]))


ONE = ((F(1), F(0)), (F(0), F(0)))

if __name__ == "__main__":
    for d in [23, 47, 71, 95, 11, 35, 59, 83]:
        u = ntype(d, 1)
        s = solutions([ONE, u, tconj(u)], True)
        print(d, d % 24, "single type n=1:", "INFEASIBLE" if not s else f"feasible {sorted(set(s))}")
