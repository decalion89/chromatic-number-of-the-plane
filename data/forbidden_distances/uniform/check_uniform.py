"""Second, independent check of the uniform certificates (certificates/interval_k*.json and certificates/tail_k*.json).

Interval certificates claim g(d, t) = z0 + z1 J0(t) + z2 J0(d t) >= 0 on [a, b] x [0, oo) with z0/(z0+z1+z2) < 1/(k-1).
Tail certificates claim h(t) = z0 + z1 J0(t) - z2 e(D t) >= 0 on [0, oo), with e(y) = min(0.4028, sqrt(2/(pi y))) and
z2 >= 0, which gives g >= h >= 0 for every d >= D (Lemma 18 of papers/forbidden-distances).

Method (different from uniform.py and tail.py, which use Taylor enclosures on t-boxes):
- the target 1/(k-1) is recomputed from k, and the intervals of a file must be contiguous and cover [a, b];
- a quick floating-point search for a counterexample, confirmed with Arb at a rational point before rejecting;
- [T, oo) by |J0(x)| <= sqrt(2/(pi x)) (Sonin-Polya);
- intervals: adaptive bisection of boxes in t AND d. On a box, x = d t lies in [d0 t0, d1 t1], and J0 on an x-interval is
  enclosed by the best of: the Arb value at the midpoint +- half-length/sqrt2 (|J0'| = |J1| <= 1/sqrt2, DLMF 10.14.1);
  the envelope sqrt(2/(pi x_lo)); the bounds -0.4028 <= J0 <= 1;
- tails: a sweep in t with steps chosen from the current value of h, the Lipschitz bound |z1|/sqrt2 for z1 J0, and
  e(D s) <= e(D t) for s >= t;
- J0 >= -0.4028 is checked here as well, by points and the Lipschitz bound on [0, 4] and the envelope beyond 4.
usage: python3 check_uniform.py FILE.json ...   (exit code 1 if any file is rejected)"""
import sys, json
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0 as fj0
from flint import arb, ctx
ctx.prec = 96
R2 = arb(2).sqrt()
PI = arb.pi()
M0 = Fr(4028, 10000)
def A(x):
    x = Fr(x); return arb(x.numerator) / arb(x.denominator)
def J(x, nu=0): return A(x).bessel_j(nu)

def j0_min_ok():
    """J0 >= -0.4028 on [0, 4] by points and Lipschitz 1/sqrt2; beyond 4 by sqrt(2/(pi x)) <= sqrt(1/(2 pi)) < 0.4028"""
    if not (arb(1) / (2 * PI)).sqrt().upper() < A(M0): return False
    t = Fr(0)
    while t < 4:
        v = J(t) + A(M0)
        if not v.lower() > 0: return False
        h = Fr(float((v.lower() * R2).mid()) * 0.9).limit_denominator(10**9)
        if h <= 0: return False
        t += h
    return True
J0MIN_OK = None

def j0_lo(xl, xh):
    """exact (arb point) lower bound of J0 on [xl, xh], 0 <= xl <= xh"""
    v = (J((xl + xh) / 2) - A((xh - xl) / 2) / R2).lower()
    c = -A(M0)
    if c > v: v = c
    if xl > 0:
        e = (-(A(2) / (PI * A(xl))).sqrt()).lower()
        if e > v: v = e
    return v

def j0_hi(xl, xh):
    v = (J((xl + xh) / 2) + A((xh - xl) / 2) / R2).upper()
    if arb(1) < v: v = arb(1)
    if xl > 0:
        e = ((A(2) / (PI * A(xl))).sqrt()).upper()
        if e < v: v = e
    return v

def scan(z, a, b, T):
    """floating-point search for a counterexample g(d, t) < 0; returns (d, t) rational, confirmed in Arb, or None"""
    zf = [float(x) for x in z]
    ts = np.arange(0.0, float(T), min(0.02, 0.2 / float(b)))
    worst = (1.0, None)
    for d in np.linspace(float(a), float(b), 9):
        g = zf[0] + zf[1] * fj0(ts) + zf[2] * fj0(d * ts)
        i = int(np.argmin(g))
        if g[i] < worst[0]: worst = (g[i], (d, ts[i]))
    if worst[0] < 0:
        d, t = worst[1]
        d = min(max(Fr(d).limit_denominator(10**6), a), b); t = Fr(t).limit_denominator(10**6)
        if (A(z[0]) + A(z[1]) * J(t) + A(z[2]) * J(d * t)).upper() < 0: return d, t
    return None

def check_interval(rec, target, max_boxes=4_000_000):
    a, b = Fr(rec['a']), Fr(rec['b'])
    z0, z1, z2 = Fr(rec['z0']), Fr(rec['z1']), Fr(rec['z2'])
    if not (0 < a < b and z0 > 0 and z0 + z1 + z2 > 0 and z0 / (z0 + z1 + z2) < target): return False, 'bound'
    T = Fr(4)
    while not (A(z0) - abs(A(z1)) * (A(2) / (PI * A(T))).sqrt()
               - abs(A(z2)) * (A(2) / (PI * A(a) * A(T))).sqrt()).lower() > 0:
        T *= 2
        if T > 2**20: return False, 'no tail'
    ce = scan((z0, z1, z2), a, b, T)
    if ce: return False, f'counterexample: g({float(ce[0])}, {float(ce[1])}) < 0'
    Z0, Z1, Z2 = A(z0), A(z1), A(z2)
    stack = [(Fr(0), T, a, b)]; boxes = 0
    while stack:
        t0, t1, d0, d1 = stack.pop(); boxes += 1
        if boxes > max_boxes: return False, f'too many boxes near t={float(t0)}'
        u = Z1 * (j0_lo(t0, t1) if z1 >= 0 else j0_hi(t0, t1))
        w = Z2 * (j0_lo(d0 * t0, d1 * t1) if z2 >= 0 else j0_hi(d0 * t0, d1 * t1))
        if (Z0 + u + w).lower() > 0: continue
        if t1 - t0 < Fr(1, 10**12): return False, f'cannot verify near t={float(t0)}, d={float(d0)}'
        if (d1 - d0) * t1 > d0 * (t1 - t0):
            dm = (d0 + d1) / 2; stack += [(t0, t1, dm, d1), (t0, t1, d0, dm)]
        else:
            tm = (t0 + t1) / 2; stack += [(tm, t1, d0, d1), (t0, tm, d0, d1)]
    return True, f'{boxes} boxes, T={float(T)}'

def e_hi(D, t):
    """an upper bound for e(D s) = min(0.4028, sqrt(2/(pi D s))) over s >= t (e is nonincreasing)"""
    if t == 0: return A(M0)
    v = (A(2) / (PI * A(D) * A(t))).sqrt().upper()
    return v if v < A(M0) else A(M0)

def check_tail(rec):
    k = rec['k']; target = Fr(1, k - 1)
    D = Fr(rec['D']); z0, z1, z2 = Fr(rec['z0']), Fr(rec['z1']), Fr(rec['z2'])
    if not (D > 0 and z0 > 0 and z2 >= 0 and z0 + z1 + z2 > 0 and z0 / (z0 + z1 + z2) < target): return False
    T = Fr(4)
    while not (A(z0) - abs(A(z1)) * (A(2) / (PI * A(T))).sqrt() - A(z2) * (A(2) / (PI * A(D) * A(T))).sqrt()).lower() > 0:
        T *= 2
        if T > 2**30: return False
    L = abs(A(z1)) / R2
    t = Fr(0); steps = 0
    while t < T:
        v = A(z0) + A(z1) * J(t) - A(z2) * e_hi(D, t)
        if not v.lower() > 0: return False
        h = Fr(float(v.lower().mid()) / float(L.upper().mid()) * 0.9).limit_denominator(10**12)
        if h <= 0: return False
        if not (v - L * A(h)).lower() > 0: return False        # z1 J0 changes by at most |z1| h/sqrt2 on [t, t+h]
        t += h; steps += 1
        if steps > 10**7: return False
    return True

if __name__ == "__main__":
    J0MIN_OK = j0_min_ok()
    allok = J0MIN_OK
    if not J0MIN_OK: print("J0 >= -0.4028 not verified")
    for f in sys.argv[1:]:
        C = json.load(open(f))
        if 'intervals' in C:
            target = Fr(1, C['k'] - 1); ivs = C['intervals']
            if 'target' in C and Fr(C['target']) != target:
                print(f"{f}: target in the file disagrees with k -> REJECTED"); allok = False; continue
            s = sorted(ivs, key=lambda r: Fr(r['a']))
            cont = all(Fr(r['a']) < Fr(r['b']) for r in ivs) and Fr(s[0]['a']) == Fr(C['a']) and \
                Fr(s[-1]['b']) == Fr(C['b']) and all(Fr(s[i]['b']) == Fr(s[i + 1]['a']) for i in range(len(s) - 1))
            bad = None
            if J0MIN_OK and cont:
                for r in s:
                    ok1, msg = check_interval(r, target)
                    if not ok1:
                        bad = r; print(f"   [{float(Fr(r['a']))}, {float(Fr(r['b']))}]: {msg}", flush=True); break
            ok = J0MIN_OK and cont and bad is None
            print(f"{f}: k={C['k']} [{float(Fr(C['a']))}, {float(Fr(C['b']))}] {len(ivs)} intervals, contiguous {cont} "
                  f"-> {'ACCEPTED' if ok else 'REJECTED'}", flush=True)
        else:
            ok = J0MIN_OK and check_tail(C)
            print(f"{f}: tail k={C['k']}, all d >= {C['D']} -> {'ACCEPTED' if ok else 'REJECTED'}", flush=True)
        allok = allok and ok
    sys.exit(0 if allok else 1)
