"""The first new components (kappa = 17/56) for N = 5, 25, 125, 625: one representative per symmetry orbit; exact
optimal point c, its decomposition c = e + z with e the nearest type point (c* or q), the values
(Re, Im)(conj(c) rho^j) mod 1 for all |j| <= k, and the '5-adic digits' at rho^{+-k}: the nearest element nu_+- of
Lambda_+- to conj(z) rho^{+-k} (mod Z[i]) and the remainder."""
from snr import *
import sys
r = F(17, 56)
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4

def nearest_lambda(w, sign):
    """nearest point of Lambda_sign to w, reduced mod Z[i] to (-1/2, 1/2]^2"""
    best = None
    g = (F(2, 5), F(sign, 5))
    for a in range(-6, 7):
        for b in range(-6, 7):
            nu = (a * g[0] - b * g[1], a * g[1] + b * g[0])
            d = (w[0] - nu[0]) ** 2 + (w[1] - nu[1]) ** 2
            if best is None or d < best[0]:
                best = (d, nu)
    nu = best[1]
    red = (nu[0] - floor(nu[0] + F(1, 2)), nu[1] - floor(nu[1] + F(1, 2)))
    return red, (w[0] - nu[0], w[1] - nu[1])

for k in range(1, KMAX + 1):
    N = 5 ** k
    L = lift_components(k, r)
    cs, qs = type_points(k)
    seen = set()
    print(f"\n== N = {N}: components with kappa 17/56 (one per symmetry orbit)")
    for P in L:
        P = clean(P)
        key = orbit_key(P, N)
        if key in seen: continue
        seen.add(key)
        kap, opt, how = component_kappa_fast(P, k)
        if kap != r: continue
        orbit = sum(1 for Q in L if orbit_key(clean(Q), N) == key)
        c = opt
        # nearest type point
        cand = [("c*", cs[0])] + [("q", q) for q in qs]
        name, e = min(cand, key=lambda t: torus_d2(c, t[1], N))
        z = (c[0] - e[0], c[1] - e[1])
        z = (z[0] - N * floor(z[0] / N + F(1, 2)), z[1] - N * floor(z[1] / N + F(1, 2)))
        vals = values(c, k)
        print(f" orbit size {orbit}: c = ({c[0]}, {c[1]}) = {name} ({e[0]}, {e[1]}) + z, z = ({z[0]}, {z[1]}) ~ ({float(z[0]):.3f}, {float(z[1]):.3f}), |z| = {math.sqrt(z[0]**2+z[1]**2):.3f}")
        print("    values:", "; ".join(f"j={j}: ({vals[j][0]}, {vals[j][1]})" for j in range(-k, k + 1)))
        for sg in (1, -1):
            A, B = rho_pow(sg * k)
            w = (A * z[0] + B * z[1], B * z[0] - A * z[1])          # conj(z) rho^{sg k}
            w = (w[0] - floor(w[0] + F(1, 2)), w[1] - floor(w[1] + F(1, 2)))
            nu, rem = nearest_lambda(w, sg)
            print(f"    at rho^{sg*k}: conj(z) rho^{sg*k} = ({w[0]}, {w[1]}) mod Z[i] = nu + rem, nu = ({nu[0]}, {nu[1]}) in Lambda_{'+' if sg>0 else '-'}, rem = ({rem[0]}, {rem[1]})")
