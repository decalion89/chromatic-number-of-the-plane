"""One certificate for all d >= D (Lemma 18 of papers/forbidden-distances).

For y > 0 put e(y) = min(mu, sqrt(2/(pi y))) with mu = 0.4028. Since J0 >= -mu on [0, oo) and |J0(x)| <= sqrt(2/(pi x)),
J0(x) >= -e(y) for every x >= y. Let
    S(t) = z0 + z1 J0(t) + sum_c y_c (2 J0(r1_c t) + J0(r3_c t)) + sum_m u_m sum_{w in M} J0(|w - q_m| t).
The first sum is over unit equilateral triangles W = {(0,0), (1,0), (1/2, sqrt3/2)} with the point p_c = centroid +
(0, oy_c), oy_c rational, so that the distances from p_c to the vertices are r1_c (twice) and r3_c:
    r1^2 = 1/3 + oy/sqrt3 + oy^2,   r3 = |1/sqrt3 - oy|.
The second is over the Moser spindle M = {0, 1, w, 1 + w, R, R w, R(1 + w)} (complex coordinates, w = e^{i pi/3},
R = (5 + i sqrt11)/6), whose unit-distance graph has independence number 2, with a rational point q_m = (qx_m, qy_m).
If z2 >= 0, y_c >= 0, u_m >= 0, s = z0 + z1 + z2 + 3 sum y_c + 7 sum u_m > 0 and
    h(t) = S(t) - z2 e(D t) >= 0   for all t > 0,
then g(d, t) = S(t) + z2 J0(d t) >= h(t) >= 0 for every d >= D and t > 0 (and g(d, 0) = s > 0), so by Proposition 17
(each triangle is a clique for {1, d}, and the spindle has at most 2 points in any set avoiding the distance 1)
m_{1,d}(R^2) <= (z0 + sum y_c + 2 sum u_m)/s for every d >= D, and chi(R^2, {1, d}) >= k for every real algebraic
d >= D when this is < 1/(k - 1).

Rigorous check of h >= 0 (Arb): on [0, T] by bisection of t-boxes [t0, t1], with J0 enclosed as in uniform.py (Taylor at
the midpoint, the envelope, -mu <= J0 <= 1) and e(D t) <= e(D t0), e being nonincreasing; the radii are enclosed by
rational intervals; on [T, oo) by z0 - |z1| sqrt(2/(pi T)) - z2 sqrt(2/(pi D T)) - (the same envelope bound for every
triangle and spindle term) > 0. The floating-point LP only proposes the coefficients.

usage: python3 tail.py k D [--cliques] [--spindles]   find and certify for all d >= D (D a fraction); writes tail_k{k}.json
       python3 tail.py check FILE.json ...            re-verify stored tail certificates (exit code 1 on failure)"""
import sys, json
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0 as fj0
from scipy.optimize import linprog
from uniform import A, lo, hi, j0_box, check_j0min, M0, PI
from flint import arb
MU = -M0                                   # 0.4028 = 1007/2500
S3 = arb(3).sqrt()
OYS = [Fr(i, 25) for i in range(-40, 41)]   # triangle positions offered to the floating-point LP
QS = [(Fr(i, 10), Fr(j, 10)) for i in range(-7, 22, 2) for j in range(-7, 22, 2)]   # spindle points offered to it

def exact(x):
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return Fr(m * 2 ** e) if e >= 0 else Fr(m, 2 ** -e)

def enclose(r):
    l, u = exact(r.lower()), exact(r.upper())
    if not 0 < l <= u: sys.exit("a radius is not separated from 0")
    return l, u

def radii(oy):
    """rational enclosures (lo, hi) of r1 and r3 for the triangle point centroid + (0, oy)"""
    o = A(oy)
    return [enclose(r) for r in ((arb(1) / 3 + o / S3 + o * o).sqrt(), abs(arb(1) / S3 - o))]

def spindle_points():
    """the seven points of the Moser spindle as Arb pairs: 0, 1, w, 1 + w, R, R w, R(1 + w)"""
    s3, s11, s33 = arb(3).sqrt(), arb(11).sqrt(), arb(33).sqrt()
    h = arb(1) / 2
    return [(arb(0), arb(0)), (arb(1), arb(0)), (h, s3 / 2), (3 * h, s3 / 2), (arb(5) / 6, s11 / 6),
            ((5 - s33) / 12, (5 * s3 + s11) / 12), ((15 - s33) / 12, (5 * s3 + 3 * s11) / 12)]

def sradii(qx, qy):
    """rational enclosures of the seven distances from the point (qx, qy) to the spindle points"""
    X, Y = A(qx), A(qy)
    return [enclose(((px - X) ** 2 + (py - Y) ** 2).sqrt()) for px, py in spindle_points()]

def _float_spindle(qx, qy):
    w = np.exp(1j * np.pi / 3); R = (5 + 1j * np.sqrt(11)) / 6
    pts = np.array([0, 1, w, 1 + w, R, R * w, R * (1 + w)])
    return np.abs(pts - complex(float(qx), float(qy)))

def float_lp(D, cliques=False, spindles=False, T=60.0):
    n = (1500, 5000, 6000) if spindles else (2000, 8000, 10000) if cliques else (3000, 16000, 30000)
    t = np.concatenate([np.linspace(1e-9, 6.0 / D, n[0]), np.linspace(6.0 / D, 8, n[1]), np.linspace(8, T, n[2])])
    e = np.minimum(float(MU), np.sqrt(2 / (np.pi * D * t)))
    E = lambda r: np.sqrt(2 / (np.pi * r * T))          # envelope row at t = T, valid for t >= T when z1 >= 0
    cols = [np.ones_like(t), fj0(t)]; env = [1, -E(1.0)]; num = [1, 0]; den = [1, 1]
    if cliques:
        for o in OYS:
            r1, r3 = np.sqrt(1/3 + float(o)/np.sqrt(3) + float(o)**2), abs(1/np.sqrt(3) - float(o))
            cols.append(2 * fj0(r1 * t) + fj0(r3 * t)); env.append(-(2 * E(r1) + E(r3))); num.append(1); den.append(3)
    if spindles:
        for qx, qy in QS:
            rs = _float_spindle(qx, qy)
            cols.append(sum(fj0(r * t) for r in rs)); env.append(-sum(E(r) for r in rs)); num.append(2); den.append(7)
    cols.append(-e); env.append(-E(D)); num.append(0); den.append(1)
    Mx = np.vstack([np.column_stack(cols), np.array(env)]); del cols
    r = linprog(num, A_ub=-Mx, b_ub=np.zeros(len(Mx)), A_eq=[den], b_eq=[1],
                bounds=[(0, None)] * len(num), method='highs')
    return r.fun, r.x

def e_upper(D, t0):
    """an upper bound for e(D t) on t >= t0"""
    if t0 == 0: return A(MU)
    s = hi((A(2) / (PI * A(D) * A(t0))).sqrt())
    return s if s < A(MU) else A(MU)

def h_lower(z, D, t0, t1, cl=(), sp=()):
    z0, z1, z2 = z
    L1, U1 = j0_box(t0, t1)
    v = A(z0) + (A(z1) * L1 if z1 >= 0 else A(z1) * U1) - A(z2) * e_upper(D, t0)
    for y, ((a1, b1), (a3, b3)) in cl:
        La, _ = j0_box(a1 * t0, b1 * t1); Lb, _ = j0_box(a3 * t0, b3 * t1)
        v = v + A(y) * (2 * La + Lb)
    for u, rs in sp:
        v = v + A(u) * sum((j0_box(a * t0, b * t1)[0] for a, b in rs), arb(0))
    return lo(v)

def tail_T(z, D, cl=(), sp=()):
    z0, z1, z2 = z
    T = Fr(10)
    env = lambda r: (A(2) / (PI * A(r) * A(T))).sqrt()
    while True:
        v = A(z0) - abs(A(z1)) * env(1) - A(z2) * env(D)
        for y, ((a1, b1), (a3, b3)) in cl:
            v = v - A(y) * (2 * env(a1) + env(a3))
        for u, rs in sp:
            v = v - A(u) * sum((env(a) for a, b in rs), arb(0))
        if lo(v) > 0: return T
        T *= 2
        if T > 10**9: return None

def verify(z, D, cl=(), sp=(), min_w=Fr(1, 10**9)):
    """(True, boxes, T) iff h >= 0 is certified on [0, oo)."""
    T = tail_T(z, D, cl, sp)
    if T is None: return False, 0, None
    stack = [(Fr(0), T)]; n = 0
    while stack:
        t0, t1 = stack.pop(); n += 1
        if h_lower(z, D, t0, t1, cl, sp) >= 0: continue
        if t1 - t0 < min_w: return False, n, T
        m = (t0 + t1) / 2; stack.append((m, t1)); stack.append((t0, m))
    return True, n, T

def bound_of(z, ys, us):
    return (z[0] + sum(ys) + 2 * sum(us)) / (sum(z) + 3 * sum(ys) + 7 * sum(us))

def admissible(k, D, z, ys=(), us=()):
    z0, z1, z2 = z
    s = z0 + z1 + z2 + 3 * sum(ys) + 7 * sum(us)
    return (D > 0 and z2 >= 0 and z0 > 0 and all(y >= 0 for y in ys) and all(u >= 0 for u in us) and s > 0
            and bound_of(z, ys, us) < Fr(1, k - 1))

def check_file(path):
    C = json.load(open(path)); k = int(C['k']); D = Fr(C['D'])
    z = (Fr(C['z0']), Fr(C['z1']), Fr(C['z2']))
    cq = [(Fr(c['oy']), Fr(c['y'])) for c in C.get('cliques', [])]
    sq = [(Fr(c['qx']), Fr(c['qy']), Fr(c['u'])) for c in C.get('spindles', [])]
    ys = [y for _, y in cq]; us = [u for _, _, u in sq]
    ok = admissible(k, D, z, ys, us) and verify(z, D, [(y, radii(oy)) for oy, y in cq],
                                                [(u, sradii(qx, qy)) for qx, qy, u in sq])[0]
    print(f"{path}: k={k}, all d >= {C['D']}, {len(cq)} triangles, {len(sq)} spindles, bound "
          f"{float(bound_of(z, ys, us)):.7f} -> {'CERTIFIED' if ok else 'FAILED'}", flush=True)
    return ok

if __name__ == "__main__":
    check_j0min()
    if sys.argv[1] == 'check':
        sys.exit(0 if all([check_file(p) for p in sys.argv[2:]]) else 1)
    k = int(sys.argv[1]); D = Fr(sys.argv[2]); target = Fr(1, k - 1)
    use_cl = '--cliques' in sys.argv[3:]; use_sp = '--spindles' in sys.argv[3:]
    v, x = float_lp(float(D), use_cl, use_sp)
    print(f"float LP at D = {float(D)}: {v:.7f} (target {float(target):.7f})", flush=True)
    if not v < float(target): sys.exit("no slack")
    nc = len(OYS) if use_cl else 0; ns = len(QS) if use_sp else 0
    z1 = Fr(x[1]).limit_denominator(10**9); z2 = Fr(x[2 + nc + ns]).limit_denominator(10**9)
    cq = [(OYS[i], Fr(x[2 + i]).limit_denominator(10**9)) for i in range(nc) if x[2 + i] > 1e-9]
    sq = [(QS[i][0], QS[i][1], Fr(x[2 + nc + i]).limit_denominator(10**9)) for i in range(ns) if x[2 + nc + i] > 1e-9]
    ys = [y for _, y in cq]; us = [u for _, _, u in sq]
    z0 = Fr(x[0]).limit_denominator(10**9) + Fr((float(target) - v) / 3).limit_denominator(10**12)
    z = (z0, z1, z2)
    if not admissible(k, D, z, ys, us): sys.exit("rounded coefficients not admissible")
    ok, n, T = verify(z, D, [(y, radii(oy)) for oy, y in cq], [(u, sradii(qx, qy)) for qx, qy, u in sq])
    bound = bound_of(z, ys, us)
    print(f"k={k} D={D}: {len(cq)} triangles, {len(sq)} spindles, bound {float(bound):.7f} < {float(target):.7f}; "
          f"{n} boxes on [0, {T}] -> {'CERTIFIED' if ok else 'FAILED'}", flush=True)
    if not ok: sys.exit(1)
    C = {'k': k, 'D': str(D), 'z0': str(z0), 'z1': str(z1), 'z2': str(z2)}
    if cq: C['cliques'] = [{'oy': str(oy), 'y': str(y)} for oy, y in cq]
    if sq: C['spindles'] = [{'qx': str(qx), 'qy': str(qy), 'u': str(u)} for qx, qy, u in sq]
    C.update({'bound': str(bound), 'T': str(T), 'boxes': n})
    json.dump(C, open(f"tail_k{k}.json", "w"))
