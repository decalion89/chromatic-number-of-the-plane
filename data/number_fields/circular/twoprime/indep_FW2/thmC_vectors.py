"""The 140 unit vectors U = G_25 V of Theorem C (d = 11, 35):
 (a) exact check that they are unit vectors, symmetric, 70 up to sign;
 (b) the 7-adic character z -> lambda(z mod w)/7, lambda(x+iy) = 2x+3y, takes values in {2,3,4,5}/7 on U
     (so kappa(U) >= 2/7; this is Proposition 3 + Theorem W+ made explicit);
 (c) a floating-point MILP (HiGHS) for kappa(U) = max_theta min_u ||theta . c_u|| over the characters of ZU
     (independent of the paper's certificates; NOT a proof).
Coordinates: z = (a + b sqrt d) + i (c + e sqrt d)  <->  (a, b, c, e) in Q^4.
"""
import sys, json, math, itertools, time
from fractions import Fraction as Fr
import numpy as np


def mul_rot(g, z):
    g1, g2 = g
    a, b, c, e = z
    return (g1 * a - g2 * c, g1 * b - g2 * e, g2 * a + g1 * c, g2 * b + g1 * e)


def build_U(d):
    def u(n):
        D = n * n + d
        return (Fr(n * n - d, D), Fr(0), Fr(0), Fr(2 * n, D))

    def conj(z):  # complex conjugation i -> -i
        a, b, c, e = z
        return (a, b, -c, -e)
    V = [(Fr(1), Fr(0), Fr(0), Fr(0))]
    for n in (1, 7, 19):
        V += [u(n), conj(u(n))]
    rho = (Fr(3, 5), Fr(4, 5))
    rots = []
    for a in range(4):
        ia = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))][a]
        for j in range(-2, 3):
            g = ia
            base = rho if j > 0 else (Fr(3, 5), Fr(-4, 5))
            for _ in range(abs(j)):
                g = (g[0] * base[0] - g[1] * base[1], g[0] * base[1] + g[1] * base[0])
            rots.append(g)
    U = []
    for g in rots:
        for v in V:
            U.append(mul_rot(g, v))
    return U


def is_unit(z, d):
    a, b, c, e = z
    return a * a + d * b * b + c * c + d * e * e == 1 and 2 * a * b + 2 * c * e == 0


def seven_adic(d, U):
    if d == 11:
        root = 2  # 2^2 = 4 = 11 mod 7
    elif d == 35:
        root = 0  # sqrt35 lies in the prime above 7 (ramified)
    vals = []
    for (a, b, c, e) in U:
        for q in (a, b, c, e):
            assert q.denominator % 7 != 0
        def red(fr):
            return fr.numerator * pow(fr.denominator, -1, 7) % 7
        x = (red(a) + root * red(b)) % 7
        y = (red(c) + root * red(e)) % 7
        vals.append((2 * x + 3 * y) % 7)
    return vals


def lattice_coords(U):
    import sympy
    from sympy.polys.matrices import DomainMatrix
    from sympy.matrices.normalforms import hermite_normal_form
    D = 1
    for z in U:
        for q in z:
            D = math.lcm(D, q.denominator)
    M = sympy.Matrix([[int(q * D) for q in z] for z in U])  # rows = vectors
    H = hermite_normal_form(M.T)  # columns span the lattice
    B = H[:, :4] if H.shape[1] >= 4 else H
    assert B.shape == (4, 4) and B.det() != 0
    # LLL on rows of B^T
    dm = DomainMatrix.from_Matrix(B.T).convert_to(sympy.ZZ)
    R = dm.lll().to_Matrix()  # rows: reduced basis
    Binv = R.T.inv()
    coords = []
    for z in U:
        col = sympy.Matrix([int(q * D) for q in z])
        cz = Binv * col
        assert all(v.is_integer for v in cz)
        coords.append([int(v) for v in cz])
    return coords, R, D


def milp_kappa(coords_half, time_limit=600):
    from scipy.optimize import milp, LinearConstraint, Bounds
    K = len(coords_half)
    # variables: theta(4), t, n_u (K)
    nv = 5 + K
    c = np.zeros(nv); c[4] = -1.0
    A = []
    lb = []
    ub = []
    for k, cu in enumerate(coords_half):
        row = np.zeros(nv); row[:4] = cu; row[4] = -1.0; row[5 + k] = -1.0
        A.append(row); lb.append(0.0); ub.append(np.inf)       # theta.c - n - t >= 0
        row = np.zeros(nv); row[:4] = cu; row[4] = 1.0; row[5 + k] = -1.0
        A.append(row); lb.append(-np.inf); ub.append(1.0)      # theta.c - n + t <= 1
    lo = np.zeros(nv); hi = np.zeros(nv)
    lo[:4] = 0; hi[:4] = 1; lo[4] = 0; hi[4] = 0.5
    for k, cu in enumerate(coords_half):
        s = sum(abs(v) for v in cu)
        lo[5 + k] = -s - 1; hi[5 + k] = s + 1
    integrality = np.zeros(nv); integrality[5:] = 1
    t0 = time.time()
    res = milp(c, constraints=LinearConstraint(np.array(A), lb, ub), integrality=integrality,
               bounds=Bounds(lo, hi), options=dict(time_limit=time_limit, disp=False, mip_rel_gap=1e-9))
    return res, time.time() - t0


if __name__ == '__main__':
    out = {}
    for d in (11, 35):
        U = build_U(d)
        assert len(U) == 140 and len(set(U)) == 140
        assert all(is_unit(z, d) for z in U)
        neg = set(tuple(-q for q in z) for z in U)
        assert neg == set(U)
        vals = seven_adic(d, U)
        rec = dict(d=d, n=len(U), unit=True, symmetric=True, seven_adic_values=sorted(set(vals)))
        print(rec, flush=True)
        coords, R, D = lattice_coords(U)
        half = []
        seen = set()
        for z, cz in zip(U, coords):
            key = tuple(cz)
            if tuple(-v for v in cz) in seen:
                continue
            seen.add(key); half.append(cz)
        rec['rank'] = 4
        rec['max_abs_coord'] = max(abs(v) for cz in half for v in cz)
        print('coords: ', len(half), 'max |coord|', rec['max_abs_coord'], flush=True)
        if len(sys.argv) > 1 and sys.argv[1] == 'milp':
            res, dt = milp_kappa(half, time_limit=float(sys.argv[2]) if len(sys.argv) > 2 else 900)
            rec['milp_status'] = res.status
            rec['milp_message'] = res.message
            rec['milp_t'] = None if res.x is None else float(res.x[4])
            rec['milp_dual_bound'] = getattr(res, 'mip_dual_bound', None)
            rec['milp_time_s'] = round(dt, 1)
            if res.x is not None:
                th = res.x[:4]
                mn = min(min((np.dot(th, cz)) % 1, 1 - (np.dot(th, cz)) % 1) for cz in half)
                rec['milp_check_min'] = float(mn)
            print(rec, flush=True)
        out[d] = rec
    json.dump(out, open('thmC_vectors_results.json', 'w'), indent=1, default=str)
