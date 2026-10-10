"""Cutting-plane LP for z0 + z1 J0(t) + z2 J0(d t) >= 0 over several d at once (float; only proposes z)."""
import numpy as np
from scipy.special import j0
from scipy.optimize import linprog
def lp_cp(ds, T=None, tol=1e-10, rounds=40):
    ds = np.atleast_1d(np.asarray(ds, float)); b = ds.max(); a = ds.min()
    T = T or max(12.0, min(40.0, 4000.0 / b))
    # initial points: fine near 0 (scale 1/b), coarse elsewhere
    t0 = np.concatenate([np.linspace(1e-7, 40.0 / b, 1500), np.linspace(40.0 / b, T, 3000)])
    P = [(d, t) for d in ds for t in t0]
    P = np.array(P)
    # fine evaluation grid for violation search
    tf = np.concatenate([np.linspace(1e-7, 60.0 / a, 6000), np.arange(60.0 / a, T, 0.05 / b)])
    for r in range(rounds):
        Mx = np.column_stack([np.ones(len(P)), j0(P[:, 1]), j0(P[:, 0] * P[:, 1])])
        res = linprog([1, 0, 0], A_ub=-Mx, b_ub=np.zeros(len(P)), A_eq=[[1, 1, 1]], b_eq=[1],
                      bounds=[(0, None), (None, None), (None, None)], method='highs')
        z = res.x
        new = []
        for d in ds:
            g = z[0] + z[1] * j0(tf) + z[2] * j0(d * tf)
            bad = np.flatnonzero(g < -tol)
            if len(bad):
                # local minima among violations
                idx = bad[(g[bad] <= np.minimum(g[np.maximum(bad - 1, 0)], g[np.minimum(bad + 1, len(g) - 1)]))]
                for i in idx[np.argsort(g[idx])[:300]]: new.append((d, tf[i]))
        if not new: return res.fun, z
        P = np.vstack([P, np.array(new)])
    return res.fun, z
if __name__ == "__main__":
    import sys, time
    for d in map(float, sys.argv[1:]):
        t = time.time(); v, z = lp_cp([d]); print(f"d={d}: LP {v:.6f} (1/LP {1/v:.4f}) z={np.round(z,5)} {time.time()-t:.1f}s", flush=True)
