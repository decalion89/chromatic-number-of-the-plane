"""Exact thresholds of Lemma C and Lemma Q when ALL constraints |j| <= k-1 are used (the note's 'r > 0.3001610'
extension), k = 2..22.  r_C(k) = 1/2 - min s,  r_Q(k) = 1/3 - min eps."""
import sys
from fractions import Fraction as Fr
from hgeom import funcs, ev, fl, gmul, gconj, rho_pow
from lpcert import exact_lp
from lp3x import lp_min
from sn import QN

def reps_pm(sg):
    g = (Fr(2, 5), Fr(sg, 5)); out = []
    for a in range(-6, 7):
        for b in range(-6, 7):
            z = gmul(g, (Fr(a), Fr(b)))
            if abs(z[0]) <= 1 and abs(z[1]) <= 1 and not (z[0].denominator == 1 and z[1].denominator == 1):
                out.append(z)
    return sorted(set(out))

def lemmaC(k, sg):
    best = None
    for nu in reps_pm(sg):
        rows = []
        for j in range(-(k - 1), k):
            for f in funcs(j):
                rows += [(f[0], f[1], Fr(-1), Fr(0)), (-f[0], -f[1], Fr(-1), Fr(0))]
        for e, f in enumerate(funcs(sg * k)):
            rows += [(f[0], f[1], Fr(-1), -nu[e]), (-f[0], -f[1], Fr(-1), nu[e])]
        res = exact_lp((0, 0, 1), rows)
        if res is not None and (best is None or res[0] < best[0]):
            best = (res[0], nu)
    return best

def lemmaQ(k, q, sg):
    best = None
    pk = gmul(gconj(q), rho_pow(sg * k))
    for nu in reps_pm(sg):
        rows = [(Fr(0), Fr(0), Fr(-1), Fr(0))]
        for j in range(-(k - 1), k):
            for f in funcs(j):
                val = ev(f, q); n = fl(val)
                rows += [(f[0], f[1], Fr(-1), n + Fr(2, 3) - val), (-f[0], -f[1], Fr(-1), val - n - Fr(1, 3))]
        for e, f in enumerate(funcs(sg * k)):
            base = (pk[e] - fl(pk[e])) + nu[e]
            rows += [(f[0], f[1], Fr(-1), Fr(2, 3) - base), (-f[0], -f[1], Fr(-1), base - Fr(1, 3))]
        res = exact_lp((0, 0, 1), rows)
        if res is not None and (best is None or res[0] < best[0]):
            best = (res[0], nu)
    return best

# cross-check the float+certificate solver against the brute-force exact solver for small k
for k in (2, 3):
    for sg in (1, -1):
        nu = reps_pm(sg)[0]
        rows = []
        for j in range(-(k - 1), k):
            for f in funcs(j):
                rows += [(f[0], f[1], Fr(-1), Fr(0)), (-f[0], -f[1], Fr(-1), Fr(0))]
        for e, f in enumerate(funcs(sg * k)):
            rows += [(f[0], f[1], Fr(-1), -nu[e]), (-f[0], -f[1], Fr(-1), nu[e])]
        a = exact_lp((0, 0, 1), rows); b = lp_min((0, 0, 1), rows)
        assert a[0] == b[0], (a[0], b[0])
print("float+certificate LP agrees with brute-force vertex enumeration on test cases")
K = int(sys.argv[1]) if len(sys.argv) > 1 else 20
for k in range(2, K + 1):
    sC = min(lemmaC(k, sg)[0] for sg in (1, -1))
    eQ = min(lemmaQ(k, q, sg)[0] for q in QN(5 ** k) for sg in (1, -1))
    rC = Fr(1, 2) - sC; rQ = Fr(1, 3) - eQ
    print(f"k={k:2d}: r_C = {float(rC):.9f} ({rC if rC.denominator < 10**12 else '...'}), r_Q = {float(rQ):.9f} ({rQ if rQ.denominator < 10**12 else '...'})")
    sys.stdout.flush()
