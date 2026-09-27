"""Rigorous check of a three-point bound from its dual solution (scripts/threepoint.py --save).

For every independent set S, its normalised frequencies v = (g, z) satisfy every constraint of the
programme and 0 <= v <= 1.  For psd multipliers Z_blk and nonnegative y, w,
    L(v) = sum_blk tr(Z_blk M_blk(v)) + sum_lin y (rhs - d.v) + sum_xi w (P0 + P.g)  >= 0,
so  |S| = F(v) <= F(v) + L(v) = 1 + const + sum_var (F_var + coef_var) v_var
          <= 1 + const + sum_var max(F_var + coef_var, 0).
With a target K0 (the programme imposes |S| >= K0), max over 0 <= v <= 1 of L(v) < 0 proves alpha < K0.

- The rotation blocks are rebuilt here, from the orbit data, in interval arithmetic (mpmath.iv, 120 bits).
- The localizing blocks have dyadic coefficients (checked), taken from build(); their corner
  n/K0 - |T| is exact.
- Each Z_blk is rounded to doubles, clipped to its nonnegative eigenvalues, shifted by 1e-9 times its
  largest eigenvalue, and proved positive definite by an exact rational LDL^T.
- The multipliers y, w of the linear constraints are re-optimised by a linear programme for the rounded
  Z_blk (any nonnegative y, w give a valid bound).

usage: python3 scripts/threepoint_verify.py certificate.npz
"""
import sys, json, time, argparse
from fractions import Fraction
import numpy as np
from mpmath import iv
from threepoint import build, character_table

iv.prec = 120


def psd_exact(Zf):
    """exact LDL^T of a symmetric matrix with float (hence rational) entries: True iff positive definite"""
    n = Zf.shape[0]
    A = [[Fraction(float(Zf[i, j])) for j in range(n)] for i in range(n)]
    for k in range(n):
        p = A[k][k]
        if p <= 0: return False
        for i in range(k + 1, n):
            f = A[i][k] / p
            if f == 0: continue
            Ai, Ak = A[i], A[k]
            for j in range(k + 1, i + 1):
                Ai[j] -= f * Ak[j]
        for i in range(k + 1, n):
            for j in range(i + 1, n):
                A[i][j] = A[j][i]
    return True


def main(path, argv):
    t0 = time.time()
    cert = np.load(path)
    meta = json.loads(str(cert['meta']))
    q, kind = meta['q'], meta['kind']
    opts = meta['options']
    K0 = meta['K0']           # |S| >= K0 tested, or None
    D = build(q, kind, verbose=False, K0=K0, tri=opts['tri'], edge=opts['edge'], trilocal=opts['trilocal'],
              pentagon=opts.get('pentagon', False))
    sol = dict(z=cert['z'], val=float(cert['value']))
    for key in ('reps', 'partner', 'shift', 'sig', 'tau', 'N', 'm', 'nrot'):
        assert D['meta'][key] == meta[key], key
    n = q * q; N = meta['N']; m = meta['m']
    reps, partner, shift = meta['reps'], meta['partner'], meta['shift']
    cls, lab, zid, unit, gcol = D['cls'], D['lab'], D['zid'], D['unit'], D['gcol']
    sig = np.array(meta['sig']); tau = np.array(meta['tau']); n0 = meta['n0']
    Qf = lambda x, y: (x * x + n0 * y * y) % q
    # sanity: sigma, tau are isometries; sigma has order N, det 1; tau det -1
    for R, dt in ((sig, 1), (tau, q - 1)):
        for (x, y) in [(1, 0), (0, 1), (1, 1), (2, 5)]:
            a, b = (R[0, 0] * x + R[0, 1] * y) % q, (R[1, 0] * x + R[1, 1] * y) % q
            assert Qf(a, b) == Qf(x, y)
        assert (R[0, 0] * R[1, 1] - R[0, 1] * R[1, 0]) % q == dt
    X, Y = np.divmod(np.arange(n), q)
    idx = lambda x, y: (x % q) * q + (y % q)
    s_img = idx((sig[0, 0] * X + sig[0, 1] * Y), (sig[1, 0] * X + sig[1, 1] * Y))
    powimg = np.zeros((N, n), int); powimg[0] = np.arange(n)
    for j in range(1, N): powimg[j] = s_img[powimg[j - 1]]
    assert np.all(s_img[powimg[N - 1]] == np.arange(n)) and len(set(powimg[:, reps[0]])) == N
    # sanity: triangle orbits have constant side-class multisets
    P = np.arange(n * n); A_, B_ = np.divmod(P, n)
    sub = lambda u, v: idx(X[u] - X[v], Y[u] - Y[v])
    key = np.sort(np.stack([cls[A_], cls[B_], cls[sub(B_, A_)]], 1), 1)
    order = np.argsort(lab, kind='stable'); lb = lab[order]; kk = key[order]
    brk = np.nonzero(np.diff(lb))[0] + 1
    starts = np.concatenate([[0], brk])
    first = kk[starts]; rep_first = np.repeat(first, np.diff(np.concatenate([starts, [len(lb)]])), axis=0)
    assert np.all(rep_first == kk), "triangle orbit with different side classes"
    del order, lb, kk, rep_first, key
    print(f"combinatorial checks passed [{time.time()-t0:.0f}s]", flush=True)

    def M(x, y, which):
        """M1(x,y) or M0(x,y) as ({var: int coef}, int const)"""
        if which == 1:
            if x == 0 and y == 0: return {}, 1
            if x == 0 or y == 0 or x == y:
                c = cls[y if x == 0 else x]
                return ({} if c == 1 else {gcol[c]: 1}), 0
            p = x * n + y
            return ({} if unit[p] else {int(zid[lab[p]]): 1}), 0
        if x == y:
            c = cls[x]
            return ({} if c == 1 else {gcol[c]: -1}), 1
        d = {}
        c = cls[sub(y, x)]
        if c != 1: d[gcol[c]] = 1
        p = x * n + y
        if not unit[p]: d[int(zid[lab[p]])] = d.get(int(zid[lab[p]]), 0) - 1
        return d, 0

    two_pi = 2 * iv.pi
    cosN = [iv.cos(two_pi * iv.mpf(t) / (2 * N)) for t in range(2 * N)]     # cos(pi t / N)
    sinN = [iv.sin(two_pi * iv.mpf(t) / (2 * N)) for t in range(2 * N)]
    r2 = iv.sqrt(2); rN = iv.sqrt(N)

    def phase(t2):   # omega^(t2/2) = exp(2 pi i t2 / (2N)), t2 integer
        t2 %= 2 * N
        return cosN[t2], sinN[t2]

    # dual vector split
    z = sol['z']
    nlin = len(D['lin']); nchar = D['Pm'].shape[0]; nv = D['nv']
    nonneg = nv + nlin + nchar
    zN = np.maximum(z[:nonneg], 0.0)
    y_lin = zN[nv:nv + nlin]; w_del = zN[nv + nlin:]
    pos = nonneg
    const = iv.mpf(0); coef = [iv.mpf(0) for _ in range(nv)]
    blk_i = 0
    shift_eps = float(argv.get('eps', 1e-9))
    for which in (1, 0):
        for k in range(N // 2 + 1):
            cols = []; seen = set()
            for a in range(m):
                if a in seen: continue
                b = partner[a]
                if b < 0: cols.append({a: (1, 0, 0)}); seen.add(a)      # (re, im, scale code: 0 -> 1, 1 -> 1/sqrt2)
                else:
                    cols.append({a: (1, 0, 1), b: (1, 0, 1)}); cols.append({a: (0, 1, 1), b: (0, -1, 1)}); seen.update((a, b))
            with_origin = (which == 1 and k == 0)
            with_empty = (which == 0 and k == 0 and K0 is not None)
            off = 1 if (with_origin or with_empty) else 0
            size = len(cols) + off
            sv = z[pos:pos + size * (size + 1) // 2]; pos += size * (size + 1) // 2
            Zf = np.zeros((size, size)); t = 0
            for j in range(size):
                for i in range(j + 1):
                    val = sv[t] if i == j else sv[t] / np.sqrt(2)
                    Zf[i, j] = Zf[j, i] = val; t += 1
            lam, V = np.linalg.eigh(Zf)
            Zp = (V * np.maximum(lam, 0)) @ V.T + shift_eps * max(1.0, lam.max()) * np.eye(size)
            Zp = (Zp + Zp.T) / 2
            assert psd_exact(Zp), ("not psd", which, k)
            # H[a][b] entries as intervals: sum_l omega^(k l + k (s_a - s_b)) M(p_a, sigma^l p_b)
            # store real and imaginary parts of the linear form
            H = {}
            for a in range(m):
                for b in range(m):
                    accr = {}; acci = {}; cr = iv.mpf(0); ci = iv.mpf(0)
                    ys = powimg[:, reps[b]]
                    d_sh2 = int(round(2 * k * (shift[a] - shift[b])))
                    for l in range(N):
                        terms, c0 = M(reps[a], int(ys[l]), which)
                        if not terms and c0 == 0: continue
                        cr_, si_ = phase(2 * k * l + d_sh2)
                        if c0: cr += c0 * cr_; ci += c0 * si_
                        for v_, c_ in terms.items():
                            accr[v_] = accr.get(v_, iv.mpf(0)) + c_ * cr_
                            acci[v_] = acci.get(v_, iv.mpf(0)) + c_ * si_
                    H[a, b] = (accr, acci, cr, ci)
            def entry(ci_, cj_):
                # sum_{a,b} conj(T_a i) T_b j H[a,b], real part (imaginary part must vanish)
                accr = {}; cr = iv.mpf(0)
                for a, (ra, ia, sa) in ci_.items():
                    for b, (rb, ib, sb) in cj_.items():
                        # w = conj(ra + i ia) (rb + i ib) = (ra rb + ia ib) + i (ra ib - ia rb)
                        wr = ra * rb + ia * ib; wi = ra * ib - ia * rb
                        sc = (iv.mpf(1) / 2) if (sa and sb) else ((1 / r2) if (sa or sb) else iv.mpf(1))
                        hr, hi, hcr, hci = H[a, b]
                        # Re(w h) = wr Re h - wi Im h
                        if wr or wi:
                            cr += sc * (wr * hcr - wi * hci)
                            for v_ in set(hr) | set(hi):
                                accr[v_] = accr.get(v_, iv.mpf(0)) + sc * (wr * hr.get(v_, 0) - wi * hi.get(v_, 0))
                return accr, cr
            for i in range(size):
                for j in range(i, size):
                    if off and (i == 0 or j == 0):
                        if i == 0 and j == 0:
                            accr, cr = {}, (iv.mpf(1) if with_origin else iv.mpf(n) / K0 - 1)
                        else:
                            cj_ = cols[j - 1]; accr = {}; cr = iv.mpf(0)
                            for b, (rb, ib, sb) in cj_.items():
                                if ib != 0: continue          # f_- has zero origin entry (same class on both orbits)
                                sc = rN * ((1 / r2) if sb else 1)
                                c = cls[reps[b]]
                                if with_origin:
                                    if c != 1: accr[gcol[c]] = accr.get(gcol[c], iv.mpf(0)) + sc * rb
                                else:
                                    cr += sc * rb
                                    if c != 1: accr[gcol[c]] = accr.get(gcol[c], iv.mpf(0)) - sc * rb
                    else:
                        accr, cr = entry(cols[i - off], cols[j - off])
                    mult = Zp[i, j] * (1 if i == j else 2)
                    if mult == 0: continue
                    const += mult * cr
                    for v_, c_ in accr.items(): coef[v_] += mult * c_
            blk_i += 1
    # extra blocks (edge / triangle localizing): exact dyadic coefficients taken from build(); the corner
    # constant 1/delta0 - k is recomputed exactly
    extra = D['blocks'][meta['nrot']:]
    for size, ent in extra:
        sv = z[pos:pos + size * (size + 1) // 2]; pos += size * (size + 1) // 2
        Zf = np.zeros((size, size)); t = 0
        for j in range(size):
            for i in range(j + 1):
                val = sv[t] if i == j else sv[t] / np.sqrt(2)
                Zf[i, j] = Zf[j, i] = val; t += 1
        lam, V = np.linalg.eigh(Zf)
        Zp = (V * np.maximum(lam, 0)) @ V.T + shift_eps * max(1.0, lam.max()) * np.eye(size)
        Zp = (Zp + Zp.T) / 2
        assert psd_exact(Zp), "extra block not psd"
        for (i, j), (acc, cst) in ent.items():
            mult = Zp[i, j] * (1 if i == j else 2)
            if mult == 0: continue
            if i == 0 and j == 0 and abs(np.real(cst) - round(np.real(cst))) > 1e-9:
                # corner k n/K0 - |T| (k = 1 for cliques, 2 for the pentagon): exact
                for kk in (1, 2, 3):
                    cc = round(kk * float(n) / K0 - np.real(cst))
                    if abs(kk * float(n) / K0 - cc - np.real(cst)) < 1e-9:
                        break
                else:
                    raise AssertionError("unrecognised corner")
                cr = kk * iv.mpf(n) / K0 - cc
            else:
                cv = np.real(cst)
                assert float(cv * 16).is_integer(), ("non-dyadic constant", cv)
                cr = iv.mpf(cv)
            const += mult * cr
            for v_, c_ in acc.items():
                c_ = np.real(c_)
                assert float(c_ * 16).is_integer(), ("non-dyadic coefficient", c_)
                coef[v_] += mult * iv.mpf(c_)
    assert pos == len(z)
    # LP polish: with the psd multipliers fixed, choose y, w >= 0 minimising the final bound
    if argv.get('polish', '1') == '1':
        import scipy.sparse as sps
        from scipy.optimize import linprog
        size_c0, Pm0, _ = character_table(q, kind)
        nlin_ = len(D['lin']); nch = Pm0.shape[0]
        Fm = np.zeros(nv)
        if not argv.get('infeasible'):
            for c in D['free']: Fm[gcol[c]] = float(D['size_c'][c])
        cz = np.array([float(iv.mpf(x).mid) if hasattr(x, 'mid') else float(x) for x in coef])
        rows_, cols_, vals_ = [], [], []
        for l, (d, rhs) in enumerate(D['lin']):
            for i, cc in d.items(): rows_.append(i); cols_.append(l); vals_.append(-float(cc))
        for j in range(nch):
            for c in D['free']:
                rows_.append(gcol[c]); cols_.append(nlin_ + j); vals_.append(float(Pm0[j, c]))
        for i in range(nv): rows_.append(i); cols_.append(nlin_ + nch + i); vals_.append(-1.0)
        Aub = sps.csr_matrix((vals_, (rows_, cols_)), shape=(nv, nlin_ + nch + nv))
        bub = -(Fm + cz)
        cobj = np.concatenate([[float(rhs) for d, rhs in D['lin']], Pm0[:, 0], np.ones(nv)])
        res = linprog(cobj, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * (nlin_ + nch + nv), method='highs')
        if res.status == 0:
            y_lin = np.maximum(res.x[:nlin_], 0); w_del = np.maximum(res.x[nlin_:nlin_ + nch], 0)
            print(f"  LP polish: predicted bound {1 + float(const.mid) + res.fun if not argv.get('infeasible') else float(const.mid) + res.fun:.6f}", flush=True)
        else:
            print(f"  LP polish failed ({res.message}); keeping the solver's multipliers", flush=True)
    # linear constraints: y (rhs - d.v) >= 0
    for yv, (d, rhs) in zip(y_lin, D['lin']):
        if yv == 0: continue
        const += iv.mpf(float(yv)) * iv.mpf(rhs)
        for v_, c_ in d.items(): coef[v_] -= iv.mpf(float(yv)) * iv.mpf(c_)
    # Delsarte: w (P0 + sum P_c g_c) >= 0, P recomputed with intervals
    size_c, Pm, _ = character_table(q, kind)
    # recompute the character table rigorously on the same orbit representatives: find a xi for each row
    Xw, Yw = X, Y
    Cl = cls
    rowsxi = []
    for r_ in range(Pm.shape[0]):
        # locate a character with this row (numerically)
        for xi in range(1, n):
            a_, b_ = divmod(xi, q)
            ph = np.cos(2 * np.pi * ((a_ * Xw + b_ * Yw) % q) / q)
            row = np.bincount(Cl, weights=ph, minlength=len(size_c))
            if np.allclose(row, Pm[r_], atol=1e-7):
                rowsxi.append(xi); break
    assert len(rowsxi) == Pm.shape[0]
    cosq = [iv.cos(two_pi * iv.mpf(t) / q) for t in range(q)]
    for wv, xi in zip(w_del, rowsxi):
        if wv == 0: continue
        a_, b_ = divmod(xi, q)
        tv = (a_ * Xw + b_ * Yw) % q
        Prow = {}
        for c_, t_ in zip(Cl, tv):
            Prow[c_] = Prow.get(c_, iv.mpf(0)) + cosq[t_]
        const += iv.mpf(float(wv)) * Prow[0]
        for c in D['free']:
            coef[gcol[c]] += iv.mpf(float(wv)) * Prow.get(c, iv.mpf(0))
    # objective F = 1 + sum |K_c| g_c
    F = [iv.mpf(0)] * nv
    for c in D['free']: F[gcol[c]] = iv.mpf(int(D['size_c'][c]))
    bound = 1 + const
    worst = iv.mpf(0)
    for i in range(nv):
        t_ = F[i] + coef[i]
        if t_.b > 0: bound += iv.mpf([0, t_.b]) if t_.a < 0 else t_; worst = max(worst.b, t_.b) if False else worst
    print(f"{kind}{q}: rigorous alpha <= {bound.b}  (float SDP value {sol['val']:.6f}) [{time.time()-t0:.0f}s]", flush=True)
    if K0 is not None:
        # the corners n/K0 - k are valid only for |S| >= K0: the bound then refutes |S| >= K0 when below K0
        print(f"  with the corners at K0 = {K0}: {'alpha < ' + str(K0) + ' PROVED' if bound.b < K0 else 'no conclusion'}", flush=True)
    # infeasibility reading (for a Farkas certificate): max over v in [0,1] of L(v)
    Lmax = const
    for i in range(nv):
        if coef[i].b > 0: Lmax += coef[i].b
    print(f"  max_v L(v) <= {Lmax.b}  ({'INFEASIBLE: no independent set of size >= ' + str(K0) if (K0 and Lmax.b < 0) else 'no infeasibility certificate'})", flush=True)
    return bound, Lmax


if __name__ == '__main__':
    path = sys.argv[1]
    argv = dict(a.split('=') for a in sys.argv[2:] if '=' in a)
    main(path, argv)
