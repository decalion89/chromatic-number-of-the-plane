"""Exact radii of the main components.
P_k^(s) = s * P_k^(1) (P_k^(1) = {x: conj(x) rho^j in [-1,1]^2, |j|<=k}); Y_k(q) for 0 < eps = 1/3 - r <= 1/30:
the cell constraints at q read  -eps <= sigma*L(y) <= 1/3 + eps  (sigma = +1 if the value at q is 1/3, -1 if 2/3);
we compute Z_k = {y : sigma*L(y) >= -1 for all constraints} (so Y_k(q) = eps*Z_k as long as the upper constraints
are inactive, which we check: max sigma*L over eps*Z_k <= 1/3 + eps) and its radius."""
from snr import *
from lemma_lp import lin, pattern

def polygon_from_halfplanes(hp, box=F(100)):
    """hp: list of (a, b, c) meaning a x + b y >= c; returns the clipped polygon (exact), starting from a big box"""
    P = [(-box, -box), (box, -box), (box, box), (-box, box)]
    for (a, b, c) in hp:
        P = clip(P, a, b, c, +1)
        if not P: return P
    return clean(P)

if __name__ == "__main__":
    for k in (1, 2, 3, 4, 5, 6):
        hp = []
        for j in range(-k, k + 1):
            for (a, b) in lin(j):
                hp.append((a, b, F(-1))); hp.append((-a, -b, F(-1)))
        P = polygon_from_halfplanes(hp)
        R2 = max(v[0] ** 2 + v[1] ** 2 for v in P)
        print(f"P_{k}^(1): {len(P)} vertices, max |x|^2 = {R2} = {float(R2):.6f}, radius {math.sqrt(R2):.6f}")

    print()
    for k in (1, 2, 3, 4, 5, 6):
        N = 5 ** k
        rads = set()
        for a0 in (1, 2):
            for b0 in (1, 2):
                q = (F(N * a0, 3), F(N * b0, 3))
                hp = []; ups = []
                for j in range(-k, k + 1):
                    p = pattern(q, j)
                    for e, (a, b) in enumerate(lin(j)):
                        sg = 1 if p[e] == F(1, 3) else -1
                        assert p[e] in (F(1, 3), F(2, 3))
                        hp.append((sg * a, sg * b, F(-1)))
                        ups.append((sg * a, sg * b))
                Z = polygon_from_halfplanes(hp)
                assert max(abs(v[0]) for v in Z) < 50, "unbounded"
                R2 = max(v[0] ** 2 + v[1] ** 2 for v in Z)
                upmax = max(a * v[0] + b * v[1] for (a, b) in ups for v in Z)   # max sigma*L over Z (eps = 1)
                # upper constraints inactive for eps*Z iff eps*upmax <= 1/3 + eps, true for eps <= 1/30 iff upmax <= 11
                assert upmax <= 11
                rads.add((R2, len(Z), upmax))
        for R2, nv, upmax in rads:
            print(f"Y_{k}(q)/eps: {nv} vertices, max |y|^2/eps^2 = {R2} = {float(R2):.6f}, radius/eps {math.sqrt(R2):.6f} (max sigma*L/eps = {upmax})")
