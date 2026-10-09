"""Exact thresholds for the unique-lift lemmas of the induction (k >= 2, N = 5^k, N' = 5^(k-1)).

(C_k) If x in P_{k-1}^(s) (i.e. conj(x) rho^j in Sq_s for |j| <= k-1) and nu in Lambda_+ \ Z[i], then
      conj(x) rho^k + nu is not in Sq_s.     [and the mirror statement with rho^{-k}, Lambda_-]
      s_C(k) := min{ s : some x, nu violate this }  (an LP in (x1, x2, s) for each nu); (C_k) holds iff s < s_C(k).
(Q_k) If q in Q_N, y in Y_{k-1}(q) (the values at rho^j, |j| <= k-1, stay in the cells of the values at q),
      nu in Lambda_+ \ Z[i], then conj(q + y) rho^k + nu is not in h + Sq_s + Z[i].   [and mirror]
      eps_Q(k) := min{ eps = 1/3 - r : violated }; (Q_k) holds iff r > 1/3 - eps_Q(k).
Every optimum is certified exactly (lp3.check_certificate).
usage: python3 lemma_lp.py KMAX
"""
import sys
from fractions import Fraction as F
from lp3 import lp3_max, check_certificate
from snr import rho_pow, floor

HALF = F(1, 2)


def lam(sign):
    """the points of Lambda_sign \\ Z[i] of modulus <= 1.3"""
    out = set()
    g = (F(2, 5), F(sign, 5))       # (2 + sign*i)/5
    for a in range(-8, 9):
        for b in range(-8, 9):
            nu = (a * g[0] - b * g[1], a * g[1] + b * g[0])
            if nu[0].denominator == 1 and nu[1].denominator == 1:
                continue
            if nu[0] ** 2 + nu[1] ** 2 <= F(169, 100):
                out.add(nu)
    return sorted(out)


def lin(j):
    """(Re, Im)(conj(y) rho^j) = (A y1 + B y2, B y1 - A y2)"""
    A, B = rho_pow(j)
    return [(A, B), (B, -A)]


def C_threshold(k, sign):
    best = None
    base = []
    for j in range(-(k - 1), k):
        for (a, b) in lin(j):
            base.append((a, b, F(-1), F(0)))
            base.append((-a, -b, F(-1), F(0)))
    for nu in lam(sign):
        rows = list(base)
        for e, (a, b) in enumerate(lin(sign * k)):
            rows.append((a, b, F(-1), -nu[e]))
            rows.append((-a, -b, F(-1), nu[e]))
        cert = lp3_max((0, 0, -1), rows)       # always feasible (x = 0, s large) and bounded (s >= 0)
        assert cert is not None
        check_certificate((0, 0, -1), rows, cert)
        s = -cert[0]
        if best is None or s < best[0]:
            best = (s, nu, cert[1])
    return best


def pattern(q, j):
    A, B = rho_pow(j)
    u = A * q[0] + B * q[1]
    v = B * q[0] - A * q[1]
    return (u - floor(u), v - floor(v))


def Q_threshold(k, sign, q):
    best = None
    base = []
    for j in range(-(k - 1), k):
        p = pattern(q, j)
        for e, (a, b) in enumerate(lin(j)):
            # 1/3 - eps <= p + L <= 2/3 + eps
            base.append((-a, -b, F(-1), p[e] - F(1, 3)))
            base.append((a, b, F(-1), F(2, 3) - p[e]))
    base.append((F(0), F(0), F(-1), F(0)))             # eps >= 0  (no upper bound: every LP is feasible)
    pk = pattern(q, sign * k)
    for nu in lam(sign):
        rows = list(base)
        for e, (a, b) in enumerate(lin(sign * k)):
            c0 = pk[e] + nu[e] - HALF
            # |L + c0| <= 1/6 + eps
            rows.append((a, b, F(-1), F(1, 6) - c0))
            rows.append((-a, -b, F(-1), F(1, 6) + c0))
        cert = lp3_max((0, 0, -1), rows)
        assert cert is not None
        check_certificate((0, 0, -1), rows, cert)
        eps = -cert[0]
        if best is None or eps < best[0]:
            best = (eps, nu, cert[1])
    return best


if __name__ == "__main__":
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for k in range(2, KMAX + 1):
        N = 5 ** k
        outC = [C_threshold(k, sg) for sg in (1, -1)]
        sC = min(o[0] for o in outC)
        print(f"k={k}: s_C = {sC} = {float(sC):.6f}, r1_C = 1/2 - s_C = {HALF - sC} = {float(HALF - sC):.6f}"
              f"   (nu = {outC[0][1]}, x* = ({outC[0][2][0]}, {outC[0][2][1]}))")
        res = []
        for a in (1, 2):
            for b in (1, 2):
                q = (F(N * a, 3), F(N * b, 3))
                for sg in (1, -1):
                    o = Q_threshold(k, sg, q)
                    if o is not None:
                        res.append((o[0], (a, b), sg, o[1], o[2]))
        eQ = min(o[0] for o in res)
        allsame = len(set(o[0] for o in res)) == 1
        print(f"      eps_Q = {eQ} = {float(eQ):.6f}, r1_Q = 1/3 - eps_Q = {F(1, 3) - eQ} = {float(F(1, 3) - eQ):.6f}"
              f"   (all 8 (q, sign) equal: {allsame}; {len(res)} feasible)")
        sys.stdout.flush()
