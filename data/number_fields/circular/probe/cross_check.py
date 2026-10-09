"""Independent floating-point cross-check of the exact enumeration: uniform random points of C/NZ[i]; every point whose
float kappa (min over the 4(2k+1) rotations of ||Re(conj(c) g)||) is >= r + 1e-9 must lie in one of the computed
components (tested exactly with the float converted to a Fraction), and every point with kappa <= r - 1e-9 must lie in
none.  Also compares the hit frequency of each component with its exact area."""
import sys, numpy as np
from snr import *
k = int(sys.argv[1]); r = F(sys.argv[2]); M = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
N = 5 ** k
comps = [clean(P) for P in lift_components(k, r)]
rng = np.random.default_rng(seed)
rho = complex(3, 4) / 5
G = np.array([e * rho ** j for j in range(-k, k + 1) for e in (1, 1j)])
def inside(P, pt):
    # exact test for a convex polygon (counterclockwise from clean()) or a point/segment
    if len(P) == 1: return P[0] == pt
    if len(P) == 2:
        return cross(P[0], P[1], pt) == 0 and min(P[0][0], P[1][0]) <= pt[0] <= max(P[0][0], P[1][0]) and min(P[0][1], P[1][1]) <= pt[1] <= max(P[0][1], P[1][1])
    return all(cross(P[i], P[(i + 1) % len(P)], pt) >= 0 for i in range(len(P)))
# bounding boxes
boxes = [(min(v[0] for v in P), max(v[0] for v in P), min(v[1] for v in P), max(v[1] for v in P)) for P in comps]
hits = [0] * len(comps); bad_in = bad_out = 0; nhigh = 0
B = 200000
for start in range(0, M, B):
    X = rng.uniform(0, N, B); Y = rng.uniform(0, N, B)
    C = X + 1j * Y
    kap = np.full(B, 0.5)
    for g in G:
        for v in ((np.conj(C) * g).real, (np.conj(C) * g).imag):
            np.minimum(kap, np.abs(v - np.round(v)), out=kap)
    rf = float(r)
    idx_hi = np.nonzero(kap >= rf + 1e-9)[0]
    idx_mid = np.nonzero(np.abs(kap - rf) < 1e-9)[0]
    nhigh += len(idx_hi)
    for i in idx_hi:
        x, y = float(X[i]), float(Y[i])
        found = False
        for t, (P, bx) in enumerate(zip(comps, boxes)):
            # translate the point by multiples of N to the box
            for dx in (-N, 0, N):
                for dy in (-N, 0, N):
                    xx, yy = x + dx, y + dy
                    if float(bx[0]) - 1e-9 <= xx <= float(bx[1]) + 1e-9 and float(bx[2]) - 1e-9 <= yy <= float(bx[3]) + 1e-9:
                        if inside(P, (F(xx), F(yy))):
                            hits[t] += 1; found = True
            if found: break
        if not found: bad_in += 1
    # points clearly outside: sample a subset near components only is costly; test those within bounding boxes
    lowmask = kap <= rf - 1e-9
    for t, (P, bx) in enumerate(zip(comps, boxes)):
        sel = np.nonzero(lowmask & (X >= float(bx[0]) % N - 1e-9) & (X <= float(bx[0]) % N + float(bx[1] - bx[0]) + 1e-9)
                         & (Y >= float(bx[2]) % N - 1e-9) & (Y <= float(bx[2]) % N + float(bx[3] - bx[2]) + 1e-9))[0]
        for i in sel:
            xx = F(float(X[i])) - (F(float(X[i])) - bx[0]) // N * N if False else F(float(X[i]))
            yy = F(float(Y[i]))
            # shift into the component's frame
            sx = bx[0] - (bx[0] % N); sy = bx[2] - (bx[2] % N)
            if inside(P, (xx + sx, yy + sy)):
                bad_out += 1
tot_area = sum(area(P) for P in comps)
print(f"N={N} r={r}: {M} samples, {nhigh} with kappa >= r+1e-9 (expected ~{float(tot_area)/N**2*M:.1f}); "
      f"not in any computed component: {bad_in}; low-kappa points inside a component: {bad_out}")
for t, P in enumerate(comps):
    exp = float(area(P)) / N ** 2 * M
    if exp > 0.5 or hits[t]:
        print(f"   component {t}: hits {hits[t]}, expected {exp:.1f}")
