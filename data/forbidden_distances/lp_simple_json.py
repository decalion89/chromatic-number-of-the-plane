"""Float LP of Oliveira Filho-Vallentin Thm 1.1 for D (no cliques), z0 >= 0, grid adapted to max(D); writes JSON for
certify_lp2.py.  usage: python3 lp_simple_json.py "1;40" out.json"""
import sys, json, numpy as np
from scipy.special import j0
from scipy.optimize import linprog
from shapes import ev
D = [ev(e) for e in sys.argv[1].split(';')]
dmax = max(D)
t = np.linspace(1e-6, 40.0, int(40 * max(1.0, dmax) * 60))
A = np.column_stack([np.ones_like(t)] + [j0(d * t) for d in D])
res = linprog(np.r_[1.0, np.zeros(len(D))], A_ub=-A, b_ub=np.zeros(len(t)), A_eq=np.ones((1, 1 + len(D))), b_eq=[1.0],
              bounds=[(0, None)] + [(None, None)] * len(D), method='highs')
print(f'D={sys.argv[1]}: LP {res.fun:.6f} (1/LP {1/res.fun:.4f}), rows {len(t)}', flush=True)
json.dump({'D': sys.argv[1], 'z': [float(v) for v in res.x], 'cliques': [], 'lp': float(res.fun)}, open(sys.argv[2], 'w'))
