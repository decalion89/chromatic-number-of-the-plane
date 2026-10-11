"""Rational certificate for m_1(R^2) <= psiT (second version: coarse program, then a polishing program).
1. The inclusion-exclusion program of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki is solved in floating point on a
   grid of t (step 0.1 up to 100, then 0.25 up to 400) (relations between pairs are left out: they repeat the pair equations).
2. Its relation multipliers w are rounded to multiples of 2^-40 and frozen; a small program in (psi0, psiT, lam, u)
   minimizes psiT subject to Phi(t) >= 1 + safety on a fine grid up to tfine, -sum lam >= 1 + tailmargin, and
   Q(S) >= 0 on every independent set.
3. psi0, lam, u are rounded to multiples of 2^-40 (if Phi dips below 1 + safety/2 on a grid up to 10000, everything is
   scaled and rounded again), and
   psiT is raised until Q(S) >= 0 on every independent set (exact integers). check_m1.py proves the result.
usage: python3 m1cert2.py CONFIG.json OUT.json [--tmax=100] [--step=0.1] [--tail=400] [--tailstep=0.25] [--tfine=1000] [--safety=2e-5] [--tailmargin=0.03]
       [--solver=ipm|clarabel] [--threads=1] [--reuse]  (--reuse: take the first program's solution from OUT.lp.npz)"""
import sys, json, time
from fractions import Fraction as F
import numpy as np
from scipy.special import j0
import highspy
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exactconf as ec, ielp
opts = dict(a.split('=', 1) for a in sys.argv[3:] if a.startswith('--') and '=' in a)
tmax = float(opts.get('--tmax', 100)); step = float(opts.get('--step', 0.1)); tfine = float(opts.get('--tfine', 1000))
tail = float(opts.get('--tail', 400)); tailstep = float(opts.get('--tailstep', 0.25))
tgrid = np.concatenate([np.arange(0, tmax + 1e-9, step), np.arange(tmax + tailstep, tail + 1e-9, tailstep)])
safety = float(opts.get('--safety', 2e-5)); tailmargin = float(opts.get('--tailmargin', 0.03))
pts, _ = ec.load(sys.argv[1]); out = sys.argv[2]
n = len(pts)
Xf = np.array([complex(float(ec.kval(re)), float(ec.kval(im))) for re, im in pts])
d2, edges = ec.graph(pts); atoms = ec.independent_sets(n, edges); pairs_all, ncl = ec.congruences(n, d2, atoms)
pairs = [(I, J) for I, J in pairs_all if len(I) >= 3]
adj = np.zeros((n, n), dtype=bool)
for i, j in edges: adj[i, j] = adj[j, i] = True
M = np.array(atoms, dtype=np.int64)
rel = [(sum(1 << a for a in I), sum(1 << b for b in J)) for I, J in pairs]
print(f"{n} points, {len(edges)} unit pairs, {len(atoms)} independent sets, {ncl} classes, {len(pairs)} relations "
      f"(sets of >= 3 points)", flush=True)
T0 = time.time()
solver = opts.get('--solver', 'ipm'); threads = int(opts.get('--threads', 1))
cache = out + '.lp.npz'                                    # the solution of the first program, for reruns
import os
if '--reuse' in sys.argv[3:] and os.path.exists(cache):
    z = np.load(cache); v, y = float(z['value']), z['y']
    assert len(y) == 2 + n + int(z['npairs']) + len(rel)
    print(f"program value {v:.8f} (read from {cache})", flush=True)
else:
    v, dt, info = ielp.solve(Xf, tmax, step, M=M, rel=rel, want_dual=True, adj=adj, solver=solver, threads=threads,
                             crossover='off', tgrid=tgrid)
    print(f"program value {v:.8f} ({info['status']}, {dt:.0f}s; grid step {step} up to {tmax}, {tailstep} up to {tail})",
          flush=True)
    y = info['y']
    np.savez(cache, value=v, y=y, npairs=sum(1 for i in range(n) for j in range(i + 1, n) if not adj[i, j]))
npair = [(i, j) for i in range(n) for j in range(i + 1, n) if not adj[i, j]]
assert len(y) == 2 + n + len(npair) + len(rel)
wc = -y[2 + n + len(npair):]
G = 2 ** 40
Wi = [int(round(float(x) * G)) for x in wc]                     # frozen relation multipliers, in units of 2^-40
# W(S) for every independent set, exact (int)
Mi = M
WS = np.zeros(len(M), dtype=object); WSi = [0] * len(M)
WSa = np.zeros(len(M), dtype=np.int64)
for (I, J), w in zip(rel, Wi):
    if w == 0: continue
    WSa[(Mi & I) == I] += w; WSa[(Mi & J) == J] -= w
assert np.abs(WSa).max() < 2 ** 62
Wf = WSa.astype(np.float64) / G
# polishing program: variables x = (psi0, psiT, lam_1..n, u_1..m)
D = np.abs(Xf[:, None] - Xf[None, :]); pr = np.array([D[i, j] for i, j in npair]); m = len(npair)
nv = 2 + n + m
bits = ((M[:, None] >> np.arange(n)[None, :]) & 1).astype(np.float64)
pi = np.array([i for i, j in npair]); pj = np.array([j for i, j in npair])
PB = bits[:, pi] * bits[:, pj]
h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', 1)
cost = np.zeros(nv); cost[1] = 1.0
lo = np.full(nv, -highspy.kHighsInf); hi = np.full(nv, highspy.kHighsInf)
h.addVars(nv, lo, hi); h.changeColsCost(nv, np.arange(nv, dtype=np.int32), cost)
def addrows(A, lower):
    A = np.asarray(A); k = A.shape[0]
    import scipy.sparse as sp
    S = sp.csr_matrix(A)
    h.addRows(k, np.asarray(lower, dtype=float), np.full(k, highspy.kHighsInf), S.nnz, S.indptr[:-1].astype(np.int32),
              S.indices.astype(np.int32), S.data)
def phirows(t):
    t = np.asarray(t, dtype=float)
    A = np.zeros((len(t), nv)); A[:, 0] = j0(t); A[:, 2:2 + n] = -1.0; A[:, 2 + n:] = -j0(t[:, None] * pr[None, :])
    addrows(A, np.full(len(t), 1 + safety))
phirows(np.arange(0, 100, 0.1))
# tail: -sum lam >= 1 + tailmargin
A = np.zeros((1, nv)); A[0, 2:2 + n] = -1.0; addrows(A, [1 + tailmargin])
# Q rows: psiT + bits.lam + PB.u >= -W(S)
for c0 in range(0, len(M), 5000):
    k = min(5000, len(M) - c0)
    A = np.zeros((k, nv)); A[:, 1] = 1.0; A[:, 2:2 + n] = bits[c0:c0 + k]; A[:, 2 + n:] = PB[c0:c0 + k]
    addrows(A, -Wf[c0:c0 + k])
tchk = np.concatenate([np.arange(0, 50, 0.005), np.arange(50, tfine + 1e-9, 0.02)])
Jchk = j0(tchk[:, None] * pr[None, :]); J0chk = j0(tchk)
h.setOptionValue('solver', 'simplex')
t1 = time.time()
for it in range(50):
    h.run()
    x = np.array(h.getSolution().col_value)
    ph = x[0] * J0chk - x[2:2 + n].sum() - Jchk @ x[2 + n:]
    bad = np.flatnonzero((ph[1:-1] < 1 + safety * 0.999) & (ph[1:-1] <= ph[:-2]) & (ph[1:-1] <= ph[2:])) + 1
    print(f"polishing round {it}: {h.modelStatusToString(h.getModelStatus())}, psiT {x[1]:.8f}, min Phi on the check grid "
          f"{ph.min():.8f}, {len(bad)} local minima below 1 + safety ({time.time() - t1:.0f}s)", flush=True)
    if len(bad) == 0 or ph.min() >= 1 + safety * 0.75: break      # the final scaling absorbs the rest
    phirows(tchk[np.unique(np.clip(np.concatenate([bad - 1, bad, bad + 1]), 0, len(tchk) - 1))])
psi0, lam, u = x[0], x[2:2 + n], x[2 + n:]
def phi(t, psi0, lam, u):
    p = psi0 * j0(t) - lam.sum()
    for r, w in zip(pr, u):
        if w != 0: p = p - w * j0(r * t)
    return p
q = lambda z: F(int(round(float(z) * G)), G)
P0 = q(psi0); L = [q(z) for z in lam]; U = [q(z) for z in u]; Wfr = [F(w, G) for w in Wi]
# Phi on a grid up to 10000 (beyond, check_m1.py uses the envelope); if it dips below 1 + safety/2, scale everything
# (Phi and Q are homogeneous) and round back onto the grid 2^-40, which moves Phi by less than 1e-9
tv = np.concatenate([np.linspace(0, 50, 2000001), np.linspace(50, 5000, 500001), np.linspace(5000, 10000, 250001)])
ph = phi(tv, float(P0), np.array([float(z) for z in L]), np.array([float(z) for z in U]))
mphi = ph.min()
print(f"min Phi on the grid {mphi:.9f} at t = {tv[ph.argmin()]:.4f}; -sum lam {-sum(float(z) for z in L):.6f}", flush=True)
if mphi < 1 + safety / 2:
    s = (1 + safety / 2) / mphi
    P0 = q(float(P0) * s); L = [q(float(z) * s) for z in L]; U = [q(float(z) * s) for z in U]
    Wfr = [q(float(w) * s) for w in Wfr]
    ph = phi(tv, float(P0), np.array([float(z) for z in L]), np.array([float(z) for z in U]))
    print(f"scaled by {s:.9f}: min Phi on the grid {ph.min():.9f}", flush=True)
den = G                                                   # every value is a multiple of 2^-40
Lq = [int(z * den) for z in L]; Uq = [int(z * den) for z in U]; Wq = [int(w * den) for w in Wfr]
if sum(abs(z) for z in Lq) + sum(abs(z) for z in Uq) + 2 * sum(abs(z) for z in Wq) >= 2 ** 62:
    raise SystemExit('the integers would overflow int64')
Qa = np.zeros(len(M), dtype=np.int64)
for i in range(n):
    if Lq[i]: Qa[(M >> i) & 1 == 1] += Lq[i]
for (i, j), w in zip(npair, Uq):
    if w: Qa[((M >> i) & 1 == 1) & ((M >> j) & 1 == 1)] += w
for (I, J), w in zip(rel, Wq):
    if w: Qa[(M & I) == I] += w; Qa[(M & J) == J] -= w
qmin = int(Qa[1:].min())
PT = F(max(0, -qmin), den) + F(1, 10 ** 9)
PT = F(-((-PT * G) // 1), G)
print(f"min over nonempty independent sets of Q - psiT = {float(F(qmin, den)):.3e}; psiT = {float(PT):.10f} "
      f"(program {v:.8f})", flush=True)
cert = dict(side='M1', points=[[ec.kstr(re), ec.kstr(im)] for re, im in pts], psi0=str(P0), psiT=str(PT),
            lam=[str(z) for z in L], u=[[i + 1, j + 1, str(w)] for (i, j), w in zip(npair, U) if w],
            wc=[[[a + 1 for a in I], [b + 1 for b in J], str(w)] for (I, J), w in zip(pairs, Wfr) if w])
json.dump(cert, open(out, 'w'))
print(f"bound {float(PT):.10f} = 1/{1 / float(PT):.5f}; wrote {out} ({time.time() - T0:.0f}s)")
