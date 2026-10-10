"""Second, independent check of the uniform certificates (certificates/interval_k*.json and certificates/tail_k*.json).

Let S(t) = z0 + z1 J0(t) + sum_c y_c (2 J0(r1_c t) + J0(r3_c t)) + sum_m u_m sum_{w in M} J0(|w - q_m| t), where the
optional clique terms are unit equilateral triangles with the point centroid + (0, oy_c): r1^2 = 1/3 + oy/sqrt3 + oy^2,
r3 = |1/sqrt3 - oy| (computed here from these formulas), and the optional spindle terms are the Moser spindle M (two
unit rhombi 0, 1, w, 1+w with w = e^{i pi/3}, the second turned by the rotation R with cos = 5/6 about 0; built here by
complex multiplication) with a rational point q_m = (qx, qy). With y_c, u_m >= 0 and s = z0 + z1 + z2 + 3 sum y_c
+ 7 sum u_m, the claimed bound is (z0 + sum y_c + 2 sum u_m)/s < 1/(k-1) (a spindle has at most 2 points in a set
avoiding distance 1).
Interval certificates claim g(d, t) = S(t) + z2 J0(d t) >= 0 on [a, b] x [0, oo).
Tail certificates claim h(t) = S(t) - z2 e(D t) >= 0 on [0, oo), with e(y) = min(0.4028, sqrt(2/(pi y))) and z2 >= 0,
which gives g >= h >= 0 for every d >= D (Lemma 18 of papers/forbidden-distances).

Method (different from uniform.py and tail.py, which use Taylor enclosures on t-boxes):
- the target 1/(k-1) is recomputed from k, and the intervals of a file must be contiguous and cover [a, b];
- a quick floating-point search for a counterexample, confirmed with Arb at a rational point before rejecting;
- [T, oo) by |J0(x)| <= sqrt(2/(pi x)) (Sonin-Polya);
- intervals: adaptive bisection of boxes in t AND d. On a box, x = d t lies in [d0 t0, d1 t1] (x = r t in
  [r_lo t0, r_hi t1] for a clique radius r in [r_lo, r_hi]), and J0 on an x-interval is enclosed by the best of: the Arb
  value at the midpoint +- half-length/sqrt2 (|J0'| = |J1| <= 1/sqrt2, DLMF 10.14.1); the envelope sqrt(2/(pi x_lo));
  the bounds -0.4028 <= J0 <= 1;
- tails: a sweep in t with steps chosen from the current value of h, the Lipschitz bound of the slow part S, and
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

S3 = arb(3).sqrt()

def q(x):
    """the exact value of an arb with radius 0 (an endpoint from .lower() or .upper()) as a Fraction"""
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return Fr(m * 2 ** e) if e >= 0 else Fr(m, 2 ** -e)

def spindle_radii(qx, qy):
    """Arb balls for the seven distances from (qx, qy) to the Moser spindle, built by complex multiplication"""
    def mul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    def add(a, b): return (a[0] + b[0], a[1] + b[1])
    one = (arb(1), arb(0)); w = (arb(1) / 2, arb(3).sqrt() / 2); R = (arb(5) / 6, arb(11).sqrt() / 6)
    pts = [(arb(0), arb(0)), one, w, add(one, w), R, mul(R, w), mul(R, add(one, w))]
    X, Y = A(qx), A(qy)
    return [((px - X) ** 2 + (py - Y) ** 2).sqrt() for px, py in pts]

def spindles_of(rec):
    """[(u, [r_1..r_7], [(lo, hi)...])] or None if a weight is negative or a distance is not separated from 0"""
    out = []
    for c in rec.get('spindles', []):
        u = Fr(c['u']); rs = spindle_radii(Fr(c['qx']), Fr(c['qy']))
        if u < 0 or not all(r.lower() > 0 for r in rs): return None
        out.append((u, rs, [(q(r.lower()), q(r.upper())) for r in rs]))
    return out

def cliques_of(rec):
    """[(y, r1, r3, r1_lo, r1_hi, r3_lo, r3_hi)] with Arb balls r1, r3 and rational enclosures; None if a weight is
    negative or a radius is not separated from 0"""
    out = []
    for c in rec.get('cliques', []):
        oy, y = A(c['oy']), Fr(c['y'])
        r1 = (arb(1) / 3 + oy / S3 + oy * oy).sqrt(); r3 = abs(arb(1) / S3 - oy)
        if y < 0 or not r3.lower() > 0: return None
        out.append((y, r1, r3, q(r1.lower()), q(r1.upper()), q(r3.lower()), q(r3.upper())))
    return out

def scan(z, a, b, T, cl=(), sp=()):
    """floating-point search for a counterexample g(d, t) < 0; returns (d, t) rational, confirmed in Arb, or None"""
    zf = [float(x) for x in z]
    ts = np.arange(0.0, float(T), min(0.02, 0.2 / float(b)))
    slow = zf[0] + zf[1] * fj0(ts)
    for y, r1, r3, *_ in cl:
        slow = slow + float(y) * (2 * fj0(float(r1.mid()) * ts) + fj0(float(r3.mid()) * ts))
    for u, rs, _ in sp:
        slow = slow + float(u) * sum(fj0(float(r.mid()) * ts) for r in rs)
    worst = (1.0, None)
    for d in np.linspace(float(a), float(b), 9):
        g = slow + zf[2] * fj0(d * ts)
        i = int(np.argmin(g))
        if g[i] < worst[0]: worst = (g[i], (d, ts[i]))
    if worst[0] < 0:
        d, t = worst[1]
        d = min(max(Fr(d).limit_denominator(10**6), a), b); t = Fr(t).limit_denominator(10**6)
        v = A(z[0]) + A(z[1]) * J(t) + A(z[2]) * J(d * t)
        for y, r1, r3, *_ in cl:
            v += A(y) * (2 * (r1 * A(t)).bessel_j(0) + (r3 * A(t)).bessel_j(0))
        for u, rs, _ in sp:
            v += A(u) * sum(((r * A(t)).bessel_j(0) for r in rs), arb(0))
        if v.upper() < 0: return d, t
    return None

def check_interval(rec, target, max_boxes=4_000_000):
    a, b = Fr(rec['a']), Fr(rec['b'])
    z0, z1, z2 = Fr(rec['z0']), Fr(rec['z1']), Fr(rec['z2'])
    cl = cliques_of(rec); sp = spindles_of(rec)
    if cl is None or sp is None: return False, 'clique or spindle'
    ysum = sum((c[0] for c in cl), Fr(0)); usum = sum((c[0] for c in sp), Fr(0))
    s = z0 + z1 + z2 + 3 * ysum + 7 * usum
    if not (0 < a < b and z0 > 0 and s > 0 and (z0 + ysum + 2 * usum) / s < target): return False, 'bound'
    def tail_ok(T):
        v = A(z0) - abs(A(z1)) * (A(2) / (PI * A(T))).sqrt() - abs(A(z2)) * (A(2) / (PI * A(a) * A(T))).sqrt()
        for y, r1, r3, *_ in cl:
            v -= A(y) * (2 * (A(2) / (PI * r1 * A(T))).sqrt() + (A(2) / (PI * r3 * A(T))).sqrt())
        for u, rs, _ in sp:
            v -= A(u) * sum(((A(2) / (PI * r * A(T))).sqrt() for r in rs), arb(0))
        return v.lower() > 0
    T = Fr(4)
    while not tail_ok(T):
        T *= 2
        if T > 2**20: return False, 'no tail'
    ce = scan((z0, z1, z2), a, b, T, cl, sp)
    if ce: return False, f'counterexample: g({float(ce[0])}, {float(ce[1])}) < 0'
    Z0, Z1, Z2 = A(z0), A(z1), A(z2)
    stack = [(Fr(0), T, a, b)]; boxes = 0
    while stack:
        t0, t1, d0, d1 = stack.pop(); boxes += 1
        if boxes > max_boxes: return False, f'too many boxes near t={float(t0)}'
        u = Z1 * (j0_lo(t0, t1) if z1 >= 0 else j0_hi(t0, t1))
        w = Z2 * (j0_lo(d0 * t0, d1 * t1) if z2 >= 0 else j0_hi(d0 * t0, d1 * t1))
        v = Z0 + u + w
        for y, _, _, a1, b1, a3, b3 in cl:
            v += A(y) * (2 * j0_lo(a1 * t0, b1 * t1) + j0_lo(a3 * t0, b3 * t1))
        for u, _, bnds in sp:
            v += A(u) * sum((j0_lo(l * t0, h * t1) for l, h in bnds), arb(0))
        if v.lower() > 0: continue
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
    """h(t) = S(t) - z2 e(D t) >= 0 for t >= 0, with S as in the module docstring (triangle radii from the formulas,
    spindle radii from the points built here)."""
    k = rec['k']; target = Fr(1, k - 1)
    D = Fr(rec['D']); z0, z1, z2 = Fr(rec['z0']), Fr(rec['z1']), Fr(rec['z2'])
    cl = []
    for c in rec.get('cliques', []):
        oy, y = A(c['oy']), Fr(c['y'])
        if y < 0: return False
        r1 = (arb(1) / 3 + oy / S3 + oy * oy).sqrt(); r3 = abs(arb(1) / S3 - oy)
        if not r3.lower() > 0: return False
        cl.append((A(y), r1, r3))
    sp = spindles_of(rec)
    if sp is None: return False
    ysum = sum(Fr(c['y']) for c in rec.get('cliques', [])); usum = sum((c[0] for c in sp), Fr(0))
    s = z0 + z1 + z2 + 3 * ysum + 7 * usum
    if not (D > 0 and z0 > 0 and z2 >= 0 and s > 0 and (z0 + ysum + 2 * usum) / s < target): return False
    def tail_ok(T):
        v = A(z0) - abs(A(z1)) * (A(2) / (PI * A(T))).sqrt() - A(z2) * (A(2) / (PI * A(D) * A(T))).sqrt()
        for y, r1, r3 in cl:
            v -= y * (2 * (A(2) / (PI * r1 * A(T))).sqrt() + (A(2) / (PI * r3 * A(T))).sqrt())
        for u, rs, _ in sp:
            v -= A(u) * sum(((A(2) / (PI * r * A(T))).sqrt() for r in rs), arb(0))
        return v.lower() > 0
    T = Fr(4)
    while not tail_ok(T):
        T *= 2
        if T > 2**30: return False
    L = (abs(A(z1)) + sum((y * (2 * r1 + r3) for y, r1, r3 in cl), arb(0))
         + sum((A(u) * sum(rs, arb(0)) for u, rs, _ in sp), arb(0))) / R2   # Lipschitz bound of the slow part
    t = Fr(0); steps = 0
    while t < T:
        tt = A(t)
        v = A(z0) + A(z1) * J(t) - A(z2) * e_hi(D, t)
        for y, r1, r3 in cl:
            v += y * (2 * (r1 * tt).bessel_j(0) + (r3 * tt).bessel_j(0))
        for u, rs, _ in sp:
            v += A(u) * sum(((r * tt).bessel_j(0) for r in rs), arb(0))
        if not v.lower() > 0: return False
        h = Fr(float(v.lower().mid()) / float(L.upper().mid()) * 0.9).limit_denominator(10**12)
        if h <= 0: return False
        if not (v - L * A(h)).lower() > 0: return False        # the slow part changes by at most L h on [t, t+h]
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
