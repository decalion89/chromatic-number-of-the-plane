"""Cutting-plane LP for g(d, t) = z0 + z1 J0(t) + sum_c y_c (2 J0(r1_c t) + J0(r3_c t)) + z2 J0(d t) >= 0 over several d at
once (float; only proposes coefficients). The optional clique terms are unit equilateral triangles with the point
centroid + (0, oy), oy in OYS (radii r1 twice and r3, see tail.py); without them g = z0 + z1 J0(t) + z2 J0(d t).
Returns (value, z) with z = (z0, z1, z2) or, with cliques, (value, z, [(oy, y), ...])."""
import numpy as np
from fractions import Fraction as Fr
from scipy.special import j0
from scipy.optimize import linprog
OYS = [Fr(i, 25) for i in range(-40, 41)]

def _clique_cols(t):
    out = []
    for o in OYS:
        fo = float(o); r1 = np.sqrt(1/3 + fo/np.sqrt(3) + fo*fo); r3 = abs(1/np.sqrt(3) - fo)
        out.append(2 * j0(r1 * t) + j0(r3 * t))
    return out

def _env(T, a, nc):
    """envelope row at t = T (valid for t >= T when z1, z2 >= 0 is not assumed: we bound |z1|, |z2| through two rows)"""
    E = lambda r: np.sqrt(2 / (np.pi * r * T))
    rows = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            rows.append([1, -s1 * E(1.0)] + [-(2 * E(np.sqrt(1/3 + float(o)/np.sqrt(3) + float(o)**2)) + E(abs(1/np.sqrt(3) - float(o)))) for o in OYS[:nc]] + [-s2 * E(a)])
    return rows

def lp_cp(ds, T=None, tol=1e-10, rounds=40, cliques=False):
    ds = np.atleast_1d(np.asarray(ds, float)); b = ds.max(); a = ds.min()
    T = T or max(12.0, min(40.0, 4000.0 / b))
    nc = len(OYS) if cliques else 0
    # initial points: fine near 0 (scale 1/b), coarse elsewhere
    t0 = np.concatenate([np.linspace(1e-7, 40.0 / b, 1500), np.linspace(40.0 / b, T, 3000)])
    P = np.array([(d, t) for d in ds for t in t0])
    # fine evaluation grid for violation search
    tf = np.concatenate([np.linspace(1e-7, 60.0 / a, 6000), np.arange(60.0 / a, T, 0.05 / b)])
    CLf = np.column_stack(_clique_cols(tf)) if nc else None
    obj = [1, 0] + [1] * nc + [0]; norm = [1, 1] + [3] * nc + [1]
    bounds = [(0, None), (None, None)] + [(0, None)] * nc + [(None, None)]
    for r in range(rounds):
        cols = [np.ones(len(P)), j0(P[:, 1])] + (_clique_cols(P[:, 1]) if nc else []) + [j0(P[:, 0] * P[:, 1])]
        Mx = np.column_stack(cols)
        if nc: Mx = np.vstack([Mx, np.array(_env(T, a, nc))])     # keep g >= 0 beyond T (clique terms can make z0 small)
        res = linprog(obj, A_ub=-Mx, b_ub=np.zeros(len(Mx)), A_eq=[norm], b_eq=[1], bounds=bounds, method='highs')
        x = res.x
        new = []
        for d in ds:
            g = x[0] + x[1] * j0(tf) + x[2 + nc] * j0(d * tf) + (CLf @ x[2:2 + nc] if nc else 0)
            bad = np.flatnonzero(g < -tol)
            if len(bad):
                idx = bad[(g[bad] <= np.minimum(g[np.maximum(bad - 1, 0)], g[np.minimum(bad + 1, len(g) - 1)]))]
                for i in idx[np.argsort(g[idx])[:300]]: new.append((d, tf[i]))
        if not new: break
        P = np.vstack([P, np.array(new)])
    z = (x[0], x[1], x[2 + nc])
    if not nc: return res.fun, np.array(z)
    return res.fun, np.array(z), [(OYS[i], x[2 + i]) for i in range(nc) if x[2 + i] > 1e-9]

if __name__ == "__main__":
    import sys, time
    cl = '--cliques' in sys.argv
    for d in [float(v) for v in sys.argv[1:] if not v.startswith('--')]:
        t = time.time(); out = lp_cp([d], cliques=cl); print(f"d={d}: LP {out[0]:.6f} (1/LP {1/out[0]:.4f}) {time.time()-t:.1f}s", flush=True)
