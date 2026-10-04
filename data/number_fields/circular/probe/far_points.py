"""r <= 3/10: the explicit points c_k = c*_N + rho^k/5 - 5^(k-1)(2 - i) of S_N^(3/10) (N = 5^k), far from all type
points.  Checks exactly: kappa(c_k) = 3/10 (min over the 4(2k+1) rotations), values at rho^{+-k}, and the distance
to c*_N + N Z[i] and to Q_N + N Z[i]."""
from snr import *
for k in range(1, 11):
    N = 5 ** k; Np = 5 ** (k - 1)
    rk = rho_pow(k)
    c = (F(N, 2) + rk[0] / 5 - 2 * Np, F(N, 2) + rk[1] / 5 + Np)
    kap = kappa_point(c, k)
    vals = values(c, k)
    cs, qs = type_points(k)
    dc = torus_d2(c, cs[0], N)
    dq = min(torus_d2(c, q, N) for q in qs)
    print(f"k={k}: kappa = {kap}; values at rho^k: {vals[k]}, at rho^-k: ({float(vals[-k][0]):.5f}, {float(vals[-k][1]):.5f}); "
          f"dist to c*-lattice = {math.sqrt(dc):.4f} (N/sqrt5 = {N/math.sqrt(5):.4f}), to Q_N-lattice = {math.sqrt(dq):.4f} "
          f"(N*sqrt2/30 = {N*math.sqrt(2)/30:.4f})")
