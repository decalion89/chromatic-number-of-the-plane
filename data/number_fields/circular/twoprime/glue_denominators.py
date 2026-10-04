"""(A) Unbounded denominators of the one-prime probe at theta = 3/10.

The bi-infinite sequence  g_j = rho^(j-t)/5 (j < t),  g_t = -(1+i)/5,  g_j = i rho^(j-t)/5 (j > t)  lies in B_{1/5}
and satisfies g_{j+1} - rho g_j in Lambda_+ (digits 0 except nu_t = -(2+i)/5, nu_{t+1} = i(2+i)/5... = (-1+2i)/5),
so it is a character of A = Z[1/5][i] with kappa = 3/10 (an 'eternal' single-glue character).  For each level K we
reconstruct exactly its point C_K mod 5^K (lifting digit by digit), check C_K in S_N^(3/10), and compute the least m
such that C_K = N eps + x with eps in (1/(6*5^m)) Z[i] and |x| <= 1/4 (as in the exact-relations lemma).  The
least m grows like K - t: no fixed finite set of types covers the probe below 10/3."""
from fractions import Fraction as Fr
import math
from gdyn import rho, rhoinv, add, sub, cmul, rho_pow, supnorm
H = (Fr(1, 2), Fr(1, 2))
def g_of(j, t):
    if j < t:  return cmul(rho_pow(j - t), (Fr(1, 5), Fr(0)))
    if j == t: return (Fr(-1, 5), Fr(-1, 5))
    return cmul(rho_pow(j - t), (Fr(0), Fr(1, 5)))
def frac(q): return q - (q.numerator // q.denominator)
def conjmul(c, z):   # conj(c) * z
    return (c[0] * z[0] + c[1] * z[1], c[0] * z[1] - c[1] * z[0])
def same_mod_Zi(p, q): return frac(p[0] - q[0]) == 0 and frac(p[1] - q[1]) == 0
def point(K, t):
    # level 0: conj(C) = h + g_0 mod Z[i]  ->  C = conj(h + g_0)
    v = add(H, g_of(0, t)); C = (v[0], -v[1])
    for k in range(1, K + 1):
        Np = 5 ** (k - 1); ok = None
        for a in range(5):
            for b in range(5):
                Cc = (C[0] + Np * a, C[1] + Np * b)
                if all(same_mod_Zi(conjmul(Cc, rho_pow(s * k)), add(H, g_of(s * k, t))) for s in (1, -1)):
                    assert ok is None; ok = Cc
        C = ok
    return C
def kappa(C, K):
    m = Fr(1, 2)
    for j in range(-K, K + 1):
        z = conjmul(C, rho_pow(j))
        for v in z: f = frac(v); m = min(m, f, 1 - f)
    return m
assert all(supnorm(g_of(j, 3)) <= Fr(1, 5) for j in range(-30, 30))
for t in (1, 2):
    for K in range(t + 1, t + 9):
        N = 5 ** K; C = point(K, t)
        assert kappa(C, K) == Fr(3, 10)
        best = None
        for m in range(0, K + 2):
            D = 6 * 5 ** m
            z = (C[0] / N, C[1] / N)
            e = (Fr(round(z[0] * D), D), Fr(round(z[1] * D), D))
            x = ((z[0] - e[0]) * N, (z[1] - e[1]) * N)
            if math.hypot(x[0], x[1]) <= 0.25:
                best = m; break
        print(f"glue at t={t}, level K={K}: C_K in S_N^(3/10) (kappa = 3/10); least m with C_K in N*(1/(6*5^m))Z[i] + B(0,1/4): {best}")
