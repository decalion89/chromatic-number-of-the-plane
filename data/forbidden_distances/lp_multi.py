"""Cutting-plane float LP of Oliveira Filho-Vallentin for several distances D (no cliques), for distances far apart:
min z0 s.t. z0 + sum z_i J0(d_i t) >= 0 for t >= 0, z0 + sum z_i = 1, z0 >= 0. The constraint is imposed on a grid,
then at the violated local minima found on a fine grid of step 0.05/max(D) up to t = 400, until none is left.
Floating point only: it chooses a certificate, which certify_lp2.py then proves.
usage: python3 lp_multi.py "1;40;1600" out.json"""
import sys, json, numpy as np
from scipy.special import j0
from scipy.optimize import linprog
from shapes import ev
D = [ev(e) for e in sys.argv[1].split(';')]; dmax = max(D)
def G(z, t):
    return z[0] + sum(z[i + 1] * j0(D[i] * t) for i in range(len(D)))
# initial grid: dense near 0, then uniform in units of 1/dmax up to 60
t = np.unique(np.concatenate([np.linspace(0, 60, 20000), np.linspace(0, 60 / dmax * 50, 20000),
                              np.geomspace(1e-6, 60, 20000)]))
for it in range(60):
    A = np.column_stack([np.ones_like(t)] + [j0(d * t) for d in D])
    res = linprog(np.r_[1.0, np.zeros(len(D))], A_ub=-A, b_ub=np.zeros(len(t)), A_eq=np.ones((1, 1 + len(D))),
                  b_eq=[1.0], bounds=[(0, None)] + [(None, None)] * len(D), method='highs')
    z = res.x
    # fine check on [0, Tmax] in chunks, step 0.05/dmax
    Tmax = 400.0
    worst = 0.0; new = []
    step = 0.05 / dmax
    for a in np.arange(0, Tmax, 2.0):
        tt = np.arange(a, a + 2.0, step)
        g = G(z, tt)
        i = np.argmin(g)
        if g[i] < -1e-10:
            # all local minima below threshold in this chunk
            m = (g[1:-1] < g[:-2]) & (g[1:-1] < g[2:]) & (g[1:-1] < -1e-10)
            idx = np.where(m)[0] + 1
            if len(idx) > 200:
                idx = idx[np.argsort(g[idx])[:200]]
            new.extend(tt[idx]); worst = min(worst, g[i])
    print(f'it {it}: LP {res.fun:.6f} (1/LP {1/res.fun:.3f}) rows {len(t)} worst {worst:.2e} new {len(new)}', flush=True)
    if not new or worst > -1e-7:
        break
    t = np.unique(np.concatenate([t, np.array(new)]))
json.dump({'D': sys.argv[1], 'z': [float(v) for v in z], 'cliques': [], 'lp': float(res.fun)}, open(sys.argv[2], 'w'))
print('z =', z)
