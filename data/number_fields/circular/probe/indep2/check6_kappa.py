# Referee check 6: kappa (= max over a component of min_gamma ||Re(conj(c) gamma)||) of the extra components,
# i.e. the least theta above which S_{5^k}^theta has only the main components (Remark (1), (3)).
# Also: the isolated points at theta = 17/56.
import sys, math, time
from fractions import Fraction as Fr
from gauss import G, RHO, dist_int
from levels import level1, lift, classify, forms, HALF
from exactlp import certified_min_s

def kappa(P, k):
    """Exact LP: min s' s.t. |l_{j,e}(c) - n_{j,e} - 1/2| <= s' for |j|<=k, with the integers of P's cell."""
    v0 = P[0]
    A = []; b = []
    for j in range(-k, k + 1):
        for (a1, a2) in forms(j):
            lv = a1 * v0[0] + a2 * v0[1]
            n = math.floor(lv)          # l(v0) in n + 1/2 + [-s,s]
            off = lv - n - HALF
            A.append([a1, a2, Fr(-1)]); b.append(-off)
            A.append([-a1, -a2, Fr(-1)]); b.append(off)
    smin, z, dual = certified_min_s(A, b)
    return HALF - smin

def scan(theta, levels_wanted, kmax):
    s = HALF - theta
    polys = level1(s)
    res = {}
    for k in range(1, kmax + 1):
        if k > 1:
            polys = lift(polys, k, s)
        if k in levels_wanted:
            ex = [P for P in polys if classify(P, k, s) == 'X']
            ks = [kappa(P, k) for P in ex]
            res[k] = (len(ex), max(ks) if ks else None, sorted(set(ks)))
            mx = res[k][1]
            print(f"theta={theta}: level {k}: {len(ex)} extra polygons, max kappa = {mx} ~ {float(mx) if mx else None}", flush=True)
            if len(set(ks)) <= 8:
                print("    distinct kappas:", [f"{x} ({float(x):.8f})" if x.denominator < 10**9 else f"{float(x):.10f}" for x in sorted(set(ks))])
    return res

if __name__ == "__main__":
    # isolated points at 17/56
    s = HALF - Fr(17, 56)
    polys = level1(s)
    for k in range(1, 5):
        if k > 1:
            polys = lift(polys, k, s)
        ex = [P for P in polys if classify(P, k, s) == 'X']
        sizes = sorted({len(P) for P in ex})
        mins = set()
        for P in ex:
            c = P[0]
            m = min(dist_int(a1 * c[0] + a2 * c[1]) for j in range(-k, k + 1) for (a1, a2) in forms(j))
            mins.add(m)
        print(f"theta=17/56 level {k}: {len(ex)} extra pieces, vertex counts {sizes}, min_gamma ||.|| at them: {mins}")
    t0 = time.time()
    scan(Fr(3005, 10000), {1, 2, 3, 4, 5, 6}, 6)        # below 3303/10981 = 0.300792
    scan(Fr(30015, 100000), {16, 17}, 17)                # below ~0.30016093
    print(f"[{time.time() - t0:.1f}s]")
