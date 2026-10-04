"""Remark (2) of Section 10: torsion characters a -> Re(conj(c) a~) of Z[1/5][i], c in (1/M)Z[i]."""
from fractions import Fraction as Fr
import math


def gmul(x, y, M):
    return ((x[0] * y[0] - x[1] * y[1]) % M, (x[0] * y[1] + x[1] * y[0]) % M)


def ginv(a, M):
    n = (a[0] * a[0] + a[1] * a[1]) % M
    ni = pow(n, -1, M)
    return ((a[0] * ni) % M, (-a[1] * ni) % M)


def group(M, gens):
    grp = {(1 % M, 0)}
    fr = [(1 % M, 0)]
    while fr:
        nf = []
        for x in fr:
            for g in gens:
                y = gmul(x, g, M)
                if y not in grp:
                    grp.add(y)
                    nf.append(y)
        fr = nf
    return grp


def kappa(M, w, grp):
    # least margin of a -> Re(conj(w/M) a~) over the group; conj(w) a = (w0 - i w1)(p + i q): Re = w0 p + w1 q
    m = M
    for (p, q) in grp:
        v = (w[0] * p + w[1] * q) % M
        m = min(m, v, M - v)
    return Fr(m, M)


for M in (41, 76):
    rho = gmul((2, 1), ginv((2, -1), M), M)
    sig = gmul((3, 2), ginv((3, -2), M), M)
    # order of rho
    x, o = rho, 1
    while x != (1 % M, 0):
        x, o = gmul(x, rho, M), o + 1
    g1 = group(M, [rho, (0, 1)])
    g2 = group(M, [rho, sig, (0, 1)])
    best = max(((kappa(M, (a, b), g1), (a, b)) for a in range(M) for b in range(M)
                if math.gcd(math.gcd(a, b), M) == 1))
    print(f"M={M}: order of rho mod M = {o}; |<i,rho>| = {len(g1)}, |<i,rho,sigma>| = {len(g2)}; "
          f"best kappa over <i,rho> among exact order M: {best[0]} at w={best[1]}")
    if M == 41:
        print("   kappa of (12+12i)/41 over <i,rho> =", kappa(41, (12, 12), g1),
              "; over <i,rho,sigma> =", kappa(41, (12, 12), g2))
    if M == 76:
        print("   best over <i,rho,sigma> for exact order 76:",
              max(kappa(M, (a, b), g2) for a in range(M) for b in range(M) if math.gcd(math.gcd(a, b), M) == 1))
