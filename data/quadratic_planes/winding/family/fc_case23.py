"""Case d = 23 (mod 24) of Theorem 4 (four-colours paper): exact decision, given Proposition 9, of whether
U = G_N u G_N u1 u G_N conj(u1) has a character into [1/3,2/3] (referee's program).
A character <=> (z, w, w') in S_N^3 with a(w + w') + b z = 0 (z = c1, w, w' determined by (c1, c2), and conversely).
With S_N = (c* + P_k + N Z[i]) u (Q_N + N Z[i]) this is: for some types, a lattice point of
shift + N Z[i] lies in lam * P_k, where lam = sum of the coefficients of the type-c variables.
usage: python3 fc_case23.py"""
from fractions import Fraction as F
import math, sys

def mul(z, w): return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])
rho = (F(3, 5), F(4, 5)); rhob = (F(3, 5), F(-4, 5)); I = (F(0), F(1)); ONE = (F(1), F(0))
def rho_pow(j):
    z = ONE
    for _ in range(abs(j)):
        z = mul(z, rho if j > 0 else rhob)
    return z

def normals(k):
    out = []
    for j in range(-k, k+1):
        r = rho_pow(j)
        out.append(r); out.append(mul(I, r))
    return out

def in_scaled_P(x, lam, nus):
    """x in lam * P_k  <=>  |<x, nu>| <= lam/6 for all nu"""
    for nu in nus:
        if abs(x[0]*nu[0] + x[1]*nu[1]) * 6 > lam:
            return False
    return True

def feasible(d, k, want_witness=False):
    N = 5**k
    a, b = (d + 1)//4, (d - 1)//2
    assert math.gcd(a, b) == 1
    nus = normals(k)
    reps = [("c", (F(N, 2), F(N, 2)))] + [("q", (F(N*al, 3), F(N*be, 3))) for al in (1, 2) for be in (1, 2)]
    R2 = F(10, 324)   # max |x|^2 on P_1 >= on P_k
    for tw, ew in reps:
        for tw2, ew2 in reps:
            for tz, ez in reps:
                lam = (a if tw == "c" else 0) + (a if tw2 == "c" else 0) + (b if tz == "c" else 0)
                sh = (a*ew[0] + a*ew2[0] + b*ez[0], a*ew[1] + a*ew2[1] + b*ez[1])
                # reduce shift mod N to [-N/2, N/2)
                s0 = sh[0] - N*math.floor(sh[0]/N + F(1, 2))
                s1 = sh[1] - N*math.floor(sh[1]/N + F(1, 2))
                if lam == 0:
                    if s0 == 0 and s1 == 0:
                        return True, (tw, tw2, tz, (0, 0))
                    continue
                rad2 = lam*lam*R2
                M = int(math.isqrt(int(rad2)) // N) + 2
                for m in range(-M, M+1):
                    for n in range(-M, M+1):
                        p = (s0 + N*m, s1 + N*n)
                        if p[0]*p[0] + p[1]*p[1] > rad2:
                            continue
                        if in_scaled_P(p, lam, nus):
                            return True, (tw, tw2, tz, p)
    return False, None

if __name__ == "__main__":
    # sharpness for N = 5, 25, 125 (d <= 1295) and N = 625 (d <= 2999), as in the Remark of Section 4
    for k, dmax in ((1, 1295), (2, 1295), (3, 1295), (4, 2999)):
        N = 5**k
        bad = []
        infeas = []
        for d in range(23, dmax+1, 24):
            f, _ = feasible(d, k)
            if not f:
                infeas.append(d)
            if (not f) != (d < F(21*N, 5)):
                bad.append(d)
        print(f"N={N}: infeasible d = {infeas[:6]}{'...' if len(infeas) > 6 else ''} (largest {max(infeas) if infeas else None});"
              f" exceptions to 'infeasible iff d < 21N/5': {bad}")
    # k = 5 and k = 6 near the threshold 21N/5
    for k in (5, 6):
        N = 5**k
        thr = F(21*N, 5)
        d0 = int(thr) - 24*30
        d0 -= (d0 - 23) % 24
        res = []
        for d in range(d0, int(thr) + 24*40, 24):
            f, wit = feasible(d, k)
            res.append((d, f))
        exc = [d for d, f in res if (not f) != (d < thr)]
        first_feasible = min(d for d, f in res if f)
        print(f"N={N}: 21N/5 = {thr}; first feasible d (d = 23 mod 24) in the window: {first_feasible};"
              f" exceptions to 'infeasible iff d < 21N/5': {exc[:5]}{'...' if len(exc) > 5 else ''} ({len(exc)} values)")
    # the binding constraint for k = 6: (N/2)(1+i)/d in P_6 needs |Im((1-i) rho^6)| * N/(2d) <= 1/6
    r6 = rho_pow(6)
    t = mul((F(1), F(-1)), r6)
    print("(1 - i) rho^6 =", t, "; (1 - i) rho =", mul((F(1), F(-1)), rho))
