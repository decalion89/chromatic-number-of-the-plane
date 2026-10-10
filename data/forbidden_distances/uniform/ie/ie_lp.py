"""Floating-point linear programs behind the inclusion-exclusion certificates (they only propose; check_ie.py proves).

A-side (two distances 1 and d >= D, Proposition 19 of the paper): maximize delta = sum kappa subject to kappa >= 0,
  C0 (f(1) = 0), the inclusion-exclusion rows of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki for X23 (IET, IE1,
  IE2, IEC over the independent atoms), the envelope row kappa(0) <= sum_{t>0} kappa(t) e(Dt) (from f(d) = 0 and
  Lemma 18), and the tangent kappa(0) >= 2 d0 delta - d0^2 (kappa(0) = delta^2); with --hit=certB.json also the
  hit-set row 1 - p_empty <= psiT_B + L_B sum_{t>0} kappa(t) e(Dt) (529 - H(t)) (Proposition 21).
B-side (--B: one distance, the hit set): maximize q = sum kappa subject to C0, IET, IE1, IE2, IEC (the program of ACMVZ).
usage: python3 ie_lp.py D d0 [--hit=certB.json] [--tmax=200] [--step=0.05] [--n0=600] --out=dual.npz
       python3 ie_lp.py --B [--tmax=200] [--step=0.05] --out=dualB.npz"""
import sys, os, json, time
from fractions import Fraction as F
import numpy as np, scipy.sparse as sp
from scipy.special import j0
from scipy.optimize import linprog
from x23 import X
opts = dict(a.split('=', 1) for a in sys.argv[1:] if a.startswith('--') and '=' in a)
BSIDE = '--B' in sys.argv
tmax = float(opts.get('--tmax', 200)); step = float(opts.get('--step', 0.05)); n0 = int(opts.get('--n0', 600))
out = opts['--out']
MU = 0.4028
n = len(X); Dm = np.abs(X[:, None] - X[None, :]); adj = np.abs(Dm - 1) < 1e-9
nb = [sum(1 << u for u in range(n) if adj[v, u]) for v in range(n)]
ind = []
def rec(start, cur, mask):
    ind.append(cur)
    for v in range(start, n):
        if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nb[v])
rec(0, 0, 0)
M = np.array(ind, dtype=np.int64); na = len(M); assert M[0] == 0
if BSIDE:
    t = np.arange(0, tmax + 1e-9, step); D = d0 = None
else:
    D = float(sys.argv[1]); d0 = float(sys.argv[2])
    tsmall = np.geomspace(1e-3 / D, 400.0 / D, n0)
    t = np.concatenate([[0.0], tsmall[tsmall < step], np.arange(step, tmax + 1e-9, step)])
nk = len(t)
pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if not adj[i, j]]
cong = [(sum(1 << (i - 1) for i in I), sum(1 << (j - 1) for j in J)) for I, J in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cong23.json')))]
rows, cols, vals, b = [], [], [], []
def add_row(kcoef=None, acols=None, avals=None, rhs=0.0):
    r = len(b)
    if kcoef is not None:
        nz = np.flatnonzero(kcoef); rows.extend([r] * len(nz)); cols.extend(nz.tolist()); vals.extend(kcoef[nz].tolist())
    if acols is not None:
        rows.extend([r] * len(acols)); cols.extend((nk + np.asarray(acols)).tolist())
        vals.extend(np.broadcast_to(avals, len(acols)).tolist())
    b.append(rhs)
add_row(kcoef=j0(t))                                                        # C0
add_row(acols=np.arange(na), avals=1.0, rhs=1.0)                            # IET
for i in range(n): add_row(kcoef=-np.ones(nk), acols=np.flatnonzero((M >> i) & 1), avals=1.0)
for i, j in pairs: add_row(kcoef=-j0(t * Dm[i, j]), acols=np.flatnonzero(((M >> i) & 1) & ((M >> j) & 1)), avals=1.0)
for I, J in cong:
    a1 = np.flatnonzero((M & I) == I); a2 = np.flatnonzero((M & J) == J)
    add_row(acols=np.concatenate([a1, a2]), avals=np.concatenate([np.ones(len(a1)), -np.ones(len(a2))]))
Aeq = sp.csr_matrix((vals, (rows, cols)), shape=(len(b), nk + na)); beq = np.array(b)
c = np.concatenate([-np.ones(nk), np.zeros(na)])
Aub = bub = None; extra = {}
if not BSIDE:
    e = np.where(t > 0, np.minimum(MU, np.sqrt(2 / (np.pi * D * np.maximum(t, 1e-300)))), 0.0)
    env = np.concatenate([-e, np.zeros(na)]); env[0] = 1.0                       # kappa(0) - sum e kappa <= 0
    tan = np.concatenate([2 * d0 * np.ones(nk), np.zeros(na)]); tan[0] -= 1.0   # 2 d0 delta - kappa(0) <= d0^2
    ub = [env, tan]; bub = [0.0, d0 * d0]
    if '--hit' in opts:
        CB = json.load(open(opts['--hit']))
        pis = (sum(F(v) for _, _, v in CB['pi']) if isinstance(CB['pi'], list)
               else len(np.flatnonzero(np.triu(adj))) * F(CB['pi']))
        psiTB = float(F(CB['psiT'])); LB = float(max(F(0), F(CB['psi0']) + pis))
        H = n + 2 * sum(j0(t * Dm[i, j]) for i in range(n) for j in range(i + 1, n))
        g = e * (n * n - H); g[0] = 0.0
        hit = np.concatenate([-LB * g, np.zeros(na)]); hit[nk] = -1.0           # -p_empty - LB sum kappa g <= psiTB - 1
        ub.append(hit); bub.append(psiTB - 1.0); extra = dict(psiTB=psiTB, LB=LB)
    Aub = sp.csr_matrix(np.vstack(ub)); bub = np.array(bub)
print(f"{'B side' if BSIDE else f'D = {D:g}, d0 = {d0}'}: {Aeq.shape[0]} equality rows, {nk} kappa + {na} atom columns",
      flush=True)
t0 = time.time()
r = linprog(c, A_eq=Aeq, b_eq=beq, A_ub=Aub, b_ub=bub, bounds=(0, None), method='highs')
assert r.status == 0, r.message
v = -r.fun
np.savez(out, t=t, yeq=r.eqlin.marginals, zub=(r.ineqlin.marginals if Aub is not None else np.zeros(0)),
         pairs=np.array(pairs), D=(D or 0.0), d0=(d0 or 0.0), value=v, pempty=r.x[nk], **extra)
print(f"value {v:.6f} = 1/{1 / v:.4f}, p_empty = {r.x[nk]:.4f}, kappa(0) = {r.x[0]:.6f}, {time.time() - t0:.0f}s; wrote {out}",
      flush=True)
