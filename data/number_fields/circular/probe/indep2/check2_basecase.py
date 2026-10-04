# Referee check 2: base case k=1 of Proposition 6.
# For each of the 25 classes lambda mod 5, the least s for which some x in B_s satisfies the conditions
# at rho^{+-1} (c = 5h + x + lambda, z = conj(x)), (a) with the centred ("minimal") shifts of the proof,
# (b) with all integer shifts. Exact vertex enumeration.
from fractions import Fraction as Fr
import itertools
from gauss import G, I, RHO, H, modZ
from exactlp import brute_min_s, solve3

def lin(fn):
    """R-linear map z -> fn(z) in C, as rows for (re, im) in variables (z1, z2)."""
    w1 = fn(G(1)); w2 = fn(G(0, 1))
    return (w1.re, w2.re), (w1.im, w2.im)

ROWS0 = lin(lambda z: z)
ROWSP = lin(lambda z: z * RHO)
ROWSM = lin(lambda z: z * RHO ** -1)

def constraints(nup, num, shp=(0, 0), shm=(0, 0)):
    A = []; b = []
    def box(rows, off):
        for (a1, a2), o in zip(rows, off):
            # |a.z + o| <= s  <=>  a.z - s <= -o  and  -a.z - s <= o
            A.append([a1, a2, Fr(-1)]); b.append(-o)
            A.append([-a1, -a2, Fr(-1)]); b.append(o)
    box(ROWS0, (Fr(0), Fr(0)))
    box(ROWSP, (nup.re + shp[0], nup.im + shp[1]))
    box(ROWSM, (num.re + shm[0], num.im + shm[1]))
    return A, b

def classify(lam):
    p = (lam / G(2, 1)).is_gint(); m = (lam / G(2, -1)).is_gint()
    return p, m

results = {}
shifts = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
for a in range(5):
    for bb in range(5):
        lam = G(a, bb)
        nup = modZ(lam.conj() * RHO); num = modZ(lam.conj() * RHO ** -1)
        smin, zmin = brute_min_s(*constraints(nup, num))
        sall = None
        for sp in shifts:
            for sm in shifts:
                v, _ = brute_min_s(*constraints(nup, num, sp, sm))
                if v is not None and (sall is None or v < sall):
                    sall = v
        results[(a, bb)] = (smin, sall, zmin)
        p, m = classify(lam)
        print(f"lambda={a}+{bb}i  nu_+ int={p!s:5} nu_- int={m!s:5}  least s (centred shifts) = {smin!s:7}  (all shifts) = {sall!s:7}  at z={[str(t) for t in zmin[:2]]}")

# orbit representatives of the proof
for lam in [(1, 0), (2, 0), (2, 2), (1, 1)]:
    print("orbit rep", lam, "least s centred =", results[lam][0], " all shifts =", results[lam][1])
# eight classes with exactly one integral nu: least s and uniqueness of the point at s = 11/56
s = Fr(11, 56)
cnt = 0
for (a, bb), (smin, sall, zmin) in results.items():
    lam = G(a, bb)
    p, m = classify(lam)
    if p != m:
        cnt += 1
        nup = modZ(lam.conj() * RHO); num = modZ(lam.conj() * RHO ** -1)
        A, b = constraints(nup, num)
        # 2D feasible set at s = 11/56: vertices
        pts = set()
        hp = [(r[0], r[1], bi - r[2] * s) for r, bi in zip(A, b)]  # a1 z1 + a2 z2 <= b + s
        for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(hp, 2):
            det = a1 * b2 - a2 * b1
            if det == 0: continue
            x = (c1 * b2 - c2 * b1) / det; y = (a1 * c2 - a2 * c1) / det
            if all(u * x + v * y <= w for u, v, w in hp):
                pts.add((x, y))
        print(f"  one-integral class {a}+{bb}i: least s = {smin}, feasible set at s=11/56 = {[(str(x), str(y)) for x, y in pts]}")
print("number of one-integral classes:", cnt)
