"""One certificate for all d >= D (Lemma 18 of papers/forbidden-distances).

For y > 0 put e(y) = min(mu, sqrt(2/(pi y))) with mu = 0.4028. Since J0 >= -mu on [0, oo) and |J0(x)| <= sqrt(2/(pi x)),
J0(x) >= -e(y) for every x >= y. Hence, if z2 >= 0, z0 + z1 + z2 > 0 and
    h(t) = z0 + z1 J0(t) - z2 e(D t) >= 0   for all t > 0,
then g(d, t) = z0 + z1 J0(t) + z2 J0(d t) >= h(t) >= 0 for every d >= D and t > 0 (and g(d, 0) = z0 + z1 + z2 > 0),
so m_{1,d}(R^2) <= z0 / (z0 + z1 + z2) for every d >= D, and chi(R^2, {1, d}) >= k for every real algebraic d >= D
when this is < 1/(k - 1).

Rigorous check of h >= 0 (Arb): on [0, T] by bisection of t-boxes [t0, t1], with J0 enclosed as in uniform.py (Taylor at
the midpoint, the envelope, -mu <= J0 <= 1) and e(D t) <= e(D t0), e being nonincreasing; on [T, oo) by
z0 - |z1| sqrt(2/(pi T)) - z2 sqrt(2/(pi D T)) > 0. The floating-point LP only proposes z.

usage: python3 tail.py k D                 find and certify z for all d >= D (D a fraction); writes tail_k{k}.json
       python3 tail.py check FILE.json ...  re-verify stored tail certificates (exit code 1 on failure)"""
import sys, json
from fractions import Fraction as Fr
import numpy as np
from scipy.special import j0 as fj0
from scipy.optimize import linprog
from uniform import A, lo, hi, j0_box, check_j0min, M0, PI
MU = -M0                                   # 0.4028 = 1007/2500

def float_lp(D, T=80.0):
    t = np.concatenate([np.linspace(1e-9, 6.0 / D, 3000), np.linspace(6.0 / D, 8, 16000), np.linspace(8, T, 30000)])
    e = np.minimum(float(MU), np.sqrt(2 / (np.pi * D * t)))
    Mx = np.column_stack([np.ones_like(t), fj0(t), -e])
    r = linprog([1, 0, 0], A_ub=-Mx, b_ub=np.zeros(len(t)), A_eq=[[1, 1, 1]], b_eq=[1],
                bounds=[(0, None), (None, None), (0, None)], method='highs')
    return r.fun, r.x

def e_upper(D, t0):
    """an upper bound for e(D t) on t >= t0"""
    if t0 == 0: return A(MU)
    s = hi((A(2) / (PI * A(D) * A(t0))).sqrt())
    return s if s < A(MU) else A(MU)

def h_lower(z, D, t0, t1):
    z0, z1, z2 = z
    L1, U1 = j0_box(t0, t1)
    return lo(A(z0) + (A(z1) * L1 if z1 >= 0 else A(z1) * U1) - A(z2) * e_upper(D, t0))

def tail_T(z, D):
    z0, z1, z2 = z
    T = Fr(10)
    while True:
        v = A(z0) - abs(A(z1)) * (A(2) / (PI * A(T))).sqrt() - A(z2) * (A(2) / (PI * A(D) * A(T))).sqrt()
        if lo(v) > 0: return T
        T *= 2
        if T > 10**9: return None

def verify(z, D, min_w=Fr(1, 10**9)):
    """(True, boxes, T) iff h >= 0 is certified on [0, oo)."""
    T = tail_T(z, D)
    if T is None: return False, 0, None
    stack = [(Fr(0), T)]; n = 0
    while stack:
        t0, t1 = stack.pop(); n += 1
        if h_lower(z, D, t0, t1) >= 0: continue
        if t1 - t0 < min_w: return False, n, T
        m = (t0 + t1) / 2; stack.append((m, t1)); stack.append((t0, m))
    return True, n, T

def admissible(k, D, z):
    z0, z1, z2 = z
    return D > 0 and z2 >= 0 and z0 > 0 and z0 + z1 + z2 > 0 and z0 / (z0 + z1 + z2) < Fr(1, k - 1)

def check_file(path):
    C = json.load(open(path)); k = int(C['k']); D = Fr(C['D'])
    z = (Fr(C['z0']), Fr(C['z1']), Fr(C['z2']))
    ok = admissible(k, D, z) and verify(z, D)[0]
    print(f"{path}: k={k}, all d >= {C['D']}, bound {float(z[0] / sum(z)):.7f} -> {'CERTIFIED' if ok else 'FAILED'}",
          flush=True)
    return ok

if __name__ == "__main__":
    check_j0min()
    if sys.argv[1] == 'check':
        sys.exit(0 if all([check_file(p) for p in sys.argv[2:]]) else 1)
    k = int(sys.argv[1]); D = Fr(sys.argv[2]); target = Fr(1, k - 1)
    v, x = float_lp(float(D))
    print(f"float LP at D = {float(D)}: {v:.7f} (target {float(target):.7f})", flush=True)
    if not v < float(target): sys.exit("no slack")
    z1 = Fr(x[1]).limit_denominator(10**9); z2 = Fr(x[2]).limit_denominator(10**9)
    delta = Fr((float(target) - v) / 3).limit_denominator(10**12)
    z0 = 1 - z1 - z2 + delta
    z = (z0, z1, z2)
    if not admissible(k, D, z): sys.exit("rounded coefficients not admissible")
    ok, n, T = verify(z, D)
    print(f"k={k} D={D}: z=({float(z0):.6f}, {float(z1):.6f}, {float(z2):.6f}) bound {float(z0 / sum(z)):.7f} "
          f"< {float(target):.7f}; {n} boxes on [0, {T}] -> {'CERTIFIED' if ok else 'FAILED'}", flush=True)
    if not ok: sys.exit(1)
    json.dump({'k': k, 'D': str(D), 'z0': str(z0), 'z1': str(z1), 'z2': str(z2), 'bound': str(z0 / sum(z)),
               'T': str(T), 'boxes': n}, open(f"tail_k{k}.json", "w"))
