"""checks: for q = (N/3)(a+bi) in Q_N, N = 5^k: conj(q) rho^j = (1/3)(-1)^k (-i)^j (a - bi) mod Z[i] for |j| <= k;
hence p_{j-2} = -p_j mod Z[i]; and conj(c*_N) rho^j in (1+i)/2 + Z[i]."""
from snr import *
ok = True
for k in range(1, 9):
    N = 5 ** k
    for a in (1, 2):
        for b in (1, 2):
            q = (F(N * a, 3), F(N * b, 3))
            for j in range(-k, k + 1):
                A, B = rho_pow(j)
                v = (A * q[0] + B * q[1], B * q[0] - A * q[1])            # conj(q) rho^j
                w = (F(a, 3), F(-b, 3))                                   # (a - b i)/3
                for _ in range(j % 4):
                    w = cmul(w, (F(0), F(-1)))                            # times (-i)
                if k % 2: w = (-w[0], -w[1])
                d = (v[0] - w[0], v[1] - w[1])
                ok &= d[0].denominator == 1 and d[1].denominator == 1
    cst = (F(N, 2), F(N, 2))
    for j in range(-k, k + 1):
        A, B = rho_pow(j)
        v = (A * cst[0] + B * cst[1], B * cst[0] - A * cst[1])
        ok &= (v[0] - F(1, 2)).denominator == 1 and (v[1] - F(1, 2)).denominator == 1
print("pattern formula and c* values verified for k <= 8:", ok)
