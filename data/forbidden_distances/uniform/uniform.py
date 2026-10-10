"""Uniform certificates for m_{1,d}(R^2) over intervals of d (the Oliveira Filho-Vallentin LP with two distances).

For an interval [a, b] and rationals z0, z1, z2 the certificate claims
    g(d, t) = z0 + z1 J0(t) + z2 J0(d t) >= 0   for all d in [a, b] and all t >= 0,
which gives m_{1,d}(R^2) <= z0 / (z0 + z1 + z2) for every d in [a, b] (Proposition 17 of papers/forbidden-distances),
hence chi(R^2, {1, d}) >= k for every real algebraic d in [a, b] when this is < 1/(k - 1).

Rigorous ingredients (all evaluations in Arb ball arithmetic, all coefficients exact fractions):
  * J0 over [m - r, m + r]: J0(m) +- (|J1(m)| r + r^2/2), since J0' = -J1 and |J0''| = |J0 - J2|/2 <= 1;
  * |J0(x)| <= sqrt(2/(pi x)) for x > 0;  -0.4028 <= J0(x) <= 1 for x >= 0 (proved by check_j0min below);
  * on [T, oo): g >= z0 - |z1| sqrt(2/(pi T)) - |z2| sqrt(2/(pi a T)) > 0.
On [0, T] the t-axis is bisected until every box [a, b] x [t0, t1] has a nonnegative lower bound for g, with J0(d t)
enclosed over the x-interval [a t0, b t1]. The floating-point cutting-plane LP (lpcp.py) only proposes z.

usage: python3 uniform.py k a b step        cover [a, b] (fractions) by certified intervals, chi >= k;
                                             writes interval_k{k}_{a}_{b}.json
       python3 uniform.py check FILE.json ...  re-verify stored interval certificates (exit code 1 on failure)"""
import sys, json, time
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0 as fj0
from scipy.optimize import linprog
from flint import arb, ctx
ctx.prec = 80
PI = arb.pi()
M0 = Fr(-4028, 10000)

def A(x):
    x = Fr(x); return arb(x.numerator) / arb(x.denominator)

def lo(x): return x.lower()
def hi(x): return x.upper()

def check_j0min():
    """J0(x) >= -0.4028 for x >= 0: interval enclosures on [0, 4], envelope sqrt(2/(pi x)) < 0.4028 on [4, oo)."""
    env4 = (A(2) / (PI * A(4))).sqrt()
    if not hi(env4) < -A(M0):
        raise SystemExit("envelope at 4 does not beat 0.4028")
    n = 4000
    for i in range(n):
        a, b = Fr(4 * i, n), Fr(4 * (i + 1), n)
        L, _ = j0_box(a, b, use_global=False)
        if not lo(L) >= A(M0):
            raise SystemExit(f"J0 >= -0.4028 not certified on [{a}, {b}]")
    return True

def j0_box(a, b, use_global=True):
    """Arb lower/upper bounds for J0 on [a, b] (Fractions, 0 <= a <= b)."""
    m = (a + b) / 2; r = (b - a) / 2
    am, ar = A(m), A(r)
    j0m = am.bessel_j(0); j1m = am.bessel_j(1)
    w = abs(j1m) * ar + ar * ar / 2
    L = lo(j0m - w); U = hi(j0m + w)
    if a > 0:
        env = hi((A(2) / (PI * A(a))).sqrt())
        if env < -L: L = -env
        if env < U: U = env
    if use_global:
        if lo(A(M0)) > L: L = lo(A(M0))
        if U > 1: U = arb(1)
    return L, U

S3 = arb(3).sqrt()

def exact(x):
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return Fr(m * 2 ** e) if e >= 0 else Fr(m, 2 ** -e)

def radii(oy):
    """rational enclosures (lo, hi) of the distances r1 (twice) and r3 from centroid + (0, oy) to the vertices of the
    unit triangle (0,0), (1,0), (1/2, sqrt3/2): r1^2 = 1/3 + oy/sqrt3 + oy^2, r3 = |1/sqrt3 - oy|"""
    o = A(oy)
    out = []
    for r in ((arb(1) / 3 + o / S3 + o * o).sqrt(), abs(arb(1) / S3 - o)):
        l, u = exact(r.lower()), exact(r.upper())
        if not 0 < l <= u: raise SystemExit("a clique radius is not separated from 0")
        out.append((l, u))
    return out

def g_lower(z, a_d, b_d, t0, t1, cl=()):
    """lower bound of z0 + z1 J0(t) + sum y (2 J0(r1 t) + J0(r3 t)) + z2 J0(d t) on [a_d, b_d] x [t0, t1]"""
    z0, z1, z2 = z
    L1, U1 = j0_box(t0, t1)
    L2, U2 = j0_box(a_d * t0, b_d * t1)
    v = A(z0) + (A(z1) * L1 if z1 >= 0 else A(z1) * U1) + (A(z2) * L2 if z2 >= 0 else A(z2) * U2)
    for y, ((a1, b1), (a3, b3)) in cl:
        La, _ = j0_box(a1 * t0, b1 * t1); Lb, _ = j0_box(a3 * t0, b3 * t1)
        v = v + A(y) * (2 * La + Lb)
    return lo(v)

def tail_T(z, a_d, cl=()):
    z0, z1, z2 = z
    T = Fr(10)
    while True:
        v = A(z0) - abs(A(z1)) * (A(2) / (PI * A(T))).sqrt() - abs(A(z2)) * (A(2) / (PI * A(a_d) * A(T))).sqrt()
        for y, ((a1, b1), (a3, b3)) in cl:
            v = v - A(y) * (2 * (A(2) / (PI * A(a1) * A(T))).sqrt() + (A(2) / (PI * A(a3) * A(T))).sqrt())
        if lo(v) > 0: return T
        T *= 2
        if T > 10**7: return None

def verify(z, a_d, b_d, cl=(), min_w=Fr(1, 10**6), stats=None):
    """True iff g >= 0 is certified on [a_d, b_d] x [0, oo)."""
    T = tail_T(z, a_d, cl)
    if T is None: return False
    stack = [(Fr(0), T)]
    nbox = 0
    while stack:
        t0, t1 = stack.pop()
        nbox += 1
        if g_lower(z, a_d, b_d, t0, t1, cl) >= 0: continue
        if t1 - t0 < min_w:
            if stats is not None: stats['fail_t'] = float(t0)
            return False
        mid = (t0 + t1) / 2
        stack.append((mid, t1)); stack.append((t0, mid))
    if stats is not None: stats['boxes'] = nbox; stats['T'] = float(T)
    return True

def bound_of(z, ys):
    return (z[0] + sum(ys)) / (sum(z) + 3 * sum(ys))

def certify_interval(a_d, b_d, target, depth=0, out=None, log=None, cliques=False, T=None):
    """Cover [a_d, b_d] by certified intervals with bound < target; returns list of records."""
    out = [] if out is None else out
    from lpcp import lp_cp
    res = lp_cp(np.linspace(float(a_d), float(b_d), 7), T=T, cliques=cliques)
    val, x = res[0], res[1]
    cq = [(oy, Fr(y).limit_denominator(10**9)) for oy, y in res[2]] if cliques else []
    slack = float(target) - val
    if slack > 1e-4:
        delta = Fr(slack / 3).limit_denominator(10**9)
        z1 = Fr(x[1]).limit_denominator(10**8); z2 = Fr(x[2]).limit_denominator(10**8)
        ys = [y for _, y in cq]
        z0 = 1 - z1 - z2 - 3 * sum(ys) + delta
        z = (z0, z1, z2)
        bound = bound_of(z, ys)
        st = {}
        if z0 > 0 and bound < target and verify(z, a_d, b_d, [(y, radii(oy)) for oy, y in cq], stats=st):
            rec = {'a': str(a_d), 'b': str(b_d), 'z0': str(z0), 'z1': str(z1), 'z2': str(z2)}
            if cq: rec['cliques'] = [{'oy': str(oy), 'y': str(y)} for oy, y in cq]
            rec.update({'bound': str(bound), 'boxes': st.get('boxes'), 'T': st.get('T')})
            out.append(rec)
            if log: print(f"  OK [{float(a_d):.5f}, {float(b_d):.5f}] LP {val:.5f} bound {float(bound):.5f} "
                          f"boxes {st.get('boxes')}", file=log, flush=True)
            return out
    if depth > 40: raise SystemExit(f"cannot certify near {float(a_d)}")
    mid = Fr((a_d + b_d) / 2).limit_denominator(10**7)
    certify_interval(a_d, mid, target, depth + 1, out, log, cliques, T)
    certify_interval(mid, b_d, target, depth + 1, out, log, cliques, T)
    return out

def check_file(path):
    """Re-verify a stored interval certificate with the method above."""
    C = json.load(open(path)); k = int(C['k']); target = Fr(1, k - 1)
    ivs = sorted(C['intervals'], key=lambda r: Fr(r['a']))
    ok = (Fr(C['target']) == target and all(Fr(r['a']) < Fr(r['b']) for r in ivs)
          and Fr(ivs[0]['a']) == Fr(C['a']) and Fr(ivs[-1]['b']) == Fr(C['b'])
          and all(Fr(ivs[i]['b']) == Fr(ivs[i + 1]['a']) for i in range(len(ivs) - 1)))
    for r in ivs:
        if not ok: break
        z = (Fr(r['z0']), Fr(r['z1']), Fr(r['z2']))
        cq = [(Fr(c['oy']), Fr(c['y'])) for c in r.get('cliques', [])]
        ys = [y for _, y in cq]
        ok = (Fr(r['a']) > 0 and z[0] > 0 and all(y >= 0 for y in ys) and sum(z) + 3 * sum(ys) > 0
              and bound_of(z, ys) < target and verify(z, Fr(r['a']), Fr(r['b']), [(y, radii(oy)) for oy, y in cq]))
    print(f"{path}: k={k} [{float(Fr(C['a']))}, {float(Fr(C['b']))}], {len(ivs)} intervals -> "
          f"{'CERTIFIED' if ok else 'FAILED'}", flush=True)
    return ok

if __name__ == "__main__":
    check_j0min()
    if sys.argv[1] == 'check':
        sys.exit(0 if all([check_file(p) for p in sys.argv[2:]]) else 1)
    k = int(sys.argv[1]); a_d = Fr(sys.argv[2]); b_d = Fr(sys.argv[3]); step = Fr(sys.argv[4])
    cliques = '--cliques' in sys.argv[5:]
    TL = [float(v[4:]) for v in sys.argv[5:] if v.startswith('--T=')]
    TL = TL[0] if TL else None
    target = Fr(1, k - 1)
    t0 = time.time(); recs = []
    x = a_d
    while x < b_d:
        y = min(x + step, b_d)
        certify_interval(x, y, target, out=recs, log=sys.stdout, cliques=cliques, T=TL)
        x = y
    json.dump({'k': k, 'target': str(target), 'a': str(a_d), 'b': str(b_d), 'intervals': recs},
              open(f"interval_k{k}_{float(a_d):g}_{float(b_d):g}.json", "w"))
    print(f"k={k}: [{float(a_d)}, {float(b_d)}] covered by {len(recs)} intervals in {time.time()-t0:.0f}s")
