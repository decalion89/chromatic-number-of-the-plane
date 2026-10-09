"""kappa at level 1 of a place with residue field F_p, p = 3 mod 4: max over linear lambda: F_{p^2} -> F_p of
min over z in mu_{p+1} of ||lambda(z)/p||.  Then the plane has a homomorphism to K_{p/q} for 1/kappa = p/q (a circular
colouring through the residue field), so chi_c <= 1/kappa."""
from fractions import Fraction as Fr
import sys
for p in [3, 7, 11, 19, 23, 31, 43]:
    # F_{p^2} = F_p(i), elements (a, b) = a + b i
    def mul(x, y): return ((x[0] * y[0] - x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)
    def pw(x, e):
        r = (1, 0)
        for _ in range(e): r = mul(r, x)
        return r
    mu = [(a, b) for a in range(p) for b in range(p) if (a, b) != (0, 0) and pw((a, b), p + 1) == (1, 0)]
    assert len(mu) == p + 1
    best = Fr(0); arg = None
    for c0 in range(p):
        for c1 in range(p):
            if (c0, c1) == (0, 0): continue
            m = min(min((c0 * a + c1 * b) % p, p - (c0 * a + c1 * b) % p) for a, b in mu)
            k = Fr(m, p)
            if k > best: best, arg = k, (c0, c1)
    print(f"p = {p}: kappa_1 = {best} ({float(best):.4f}), chi_c <= {1/best if best else 'inf'} = {float(1/best) if best else 0:.4f}, lambda = {arg}")
