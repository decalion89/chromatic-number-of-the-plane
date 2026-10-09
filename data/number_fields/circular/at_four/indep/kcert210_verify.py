"""Independent verification of the certificate for Remark (2): the character
xi(u) = th . (a, b, c, e) mod 1 on (1/210)Z[sqrt59]^2, u = ((a + b sqrt59)/210, (c + e sqrt59)/210),
restricted to Z U_210, keeps every unit vector at distance >= 15/59 from Z."""
import sys
from fractions import Fraction as Fr
from math import isqrt
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "quadratic_planes", "winding"))
from kappaD import units  # the repository's enumerator (read only), for a cross-check of the list
d, D = 59, 210
B = [[0, 0, 2, 0], [-2, 0, 0, 0], [0, -1, 1, 0], [-1, 0, 0, 1]]   # basis of Z U (rows)
th = [Fr(29, 59), Fr(30, 59), Fr(44, 59), Fr(15, 59)]          # values of xi on the basis
# solve B^T-system: find thhat in Q^4 with B[i] . thhat = th[i]
import itertools
A = [[Fr(x) for x in row] + [th[i]] for i, row in enumerate(B)]
for i in range(4):
    p = next(k for k in range(i, 4) if A[k][i] != 0); A[i], A[p] = A[p], A[i]
    for k in range(4):
        if k != i and A[k][i] != 0:
            f = A[k][i] / A[i][i]; A[k] = [x - f * y for x, y in zip(A[k], A[i])]
thhat = [A[i][4] / A[i][i] for i in range(4)]
print("character on Z^4 coordinates (a,b,c,e):", thhat)
# my own enumeration of the unit vectors
U = set()
for b in range(-isqrt(D * D // d), isqrt(D * D // d) + 1):
    for e in range(-isqrt(D * D // d), isqrt(D * D // d) + 1):
        R = D * D - d * (b * b + e * e)
        if R < 0: continue
        for a in range(-D, D + 1):
            c2 = R - a * a
            if c2 < 0: continue
            c = isqrt(c2)
            if c * c == c2:
                for cc in (c, -c):
                    if a * b + cc * e == 0: U.add((a, b, cc, e))
assert U == set(units(d, D)), "unit lists differ"
print(len(U), "unit vectors (same list as kappaD.units)")
# every u is in Z U (integer coordinates in B) -- check, and compute the margin
mm = Fr(1)
for u in U:
    v = sum(t * x for t, x in zip(thhat, u)); v -= v.numerator // v.denominator
    mm = min(mm, v, 1 - v)
print("least margin over U_210:", mm, ">", Fr(1, 4), ":", mm > Fr(1, 4))
# B generates Z U: each basis row is an integer combination of U? check via index: |det B| = index of Z U in Z^4
def det(M):
    M = [[Fr(x) for x in r] for r in M]; n = len(M); dt = Fr(1)
    for i in range(n):
        p = next((k for k in range(i, n) if M[k][i] != 0), None)
        if p is None: return 0
        if p != i: M[i], M[p] = M[p], M[i]; dt = -dt
        dt *= M[i][i]
        for k in range(i + 1, n):
            f = M[k][i] / M[i][i]; M[k] = [x - f * y for x, y in zip(M[k], M[i])]
    return dt
print("det B =", det(B))
# gcd of all 4x4 minors of the U-matrix = index of Z U in Z^4 (Smith form); sample-free: full computation
from math import gcd
Ul = sorted(U); g = 0
for comb in itertools.combinations(range(len(Ul)), 4):
    g = gcd(g, int(det([Ul[i] for i in comb])))
    if g == abs(int(det(B))): break
print("gcd of 4x4 minors of U (stops when it reaches |det B|):", g)
