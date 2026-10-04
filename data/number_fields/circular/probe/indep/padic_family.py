"""5-adic family c_k = c*_N + rho^k/5 - 5^(k-1)(2 - i), N = 5^k: membership in S_N^(3/10), kappa, distance to type
points.  Also the extra points at r = 17/56 (N = 5, 25): denominators of points and of their values."""
from fractions import Fraction as Fr
from math import sqrt
from hgeom import rho_pow, gmul, gconj, fl
from sn import GN_brute, GN_rho, kappa, cstar, QN, Pk, Yk, lift_levels, classify, make_refs

def d2_torus(p, T, N):
    best = None
    for m0 in range(-2, 3):
        for m1 in range(-2, 3):
            dx = p[0] - T[0] - N * (round((p[0] - T[0]) / N) + m0)
            dy = p[1] - T[1] - N * (round((p[1] - T[1]) / N) + m1)
            d = dx * dx + dy * dy
            if best is None or d < best:
                best = d
    return best

for k in range(1, 11):
    N = 5 ** k
    rk = rho_pow(k)
    c = (Fr(N, 2) + rk[0] / 5 - 5 ** (k - 1) * 2, Fr(N, 2) + rk[1] / 5 + 5 ** (k - 1))
    G = GN_brute(N) if k <= 5 else GN_rho(k)
    kap = kappa(c, G)
    dc = d2_torus(c, cstar(N), N)
    dq = min(d2_torus(c, q, N) for q in QN(N))
    # inside the main components at r = 3/10?  (s = 1/5, eps = 1/30): test membership of c - c*_N, c - q in the polygons
    print(f"k={k}: kappa(c_k) = {kap}; in S_N^(3/10): {kap >= Fr(3,10)}; dist to c*_N+NZ[i] = {sqrt(dc):.4f}, to Q_N+NZ[i] = {sqrt(dq):.4f}")

print("radius of P_k^(1/5) <= (1/5)sqrt(10)/3 =", (1/5)*sqrt(10)/3, "; of Y_k at eps = 1/30 <= sqrt(2)/30 =", sqrt(2)/30)

r = Fr(17, 56)
lev = lift_levels(2, r)
for k in (1, 2):
    N = 5 ** k
    refs = make_refs(k, r)
    G = GN_brute(N)
    for P in lev[k]:
        if classify(P, k, r, refs) is None:
            v = P.V[0]
            vals = []
            for g in G:
                t = v[0] * g[0] + v[1] * g[1]
                vals.append(t - fl(t))
            dens = sorted({x.denominator for x in vals})
            mind = min(min(x, 1 - x) for x in vals)
            print(f"N={N} extra point ({v[0]}, {v[1]}): point denominators {v[0].denominator},{v[1].denominator}; value denominators {dens}; min ||.|| = {mind}; #values equal to 17/56 or 39/56: {sum(1 for x in vals if x in (Fr(17,56), Fr(39,56)))}")
