"""Rational certificate for m_1(R^3) <= B: the inclusion-exclusion program in R^3 (kernel sin(s)/s) with the atom
kappa(0) = delta^2, for an exact configuration (conf3.py format).
Certificate: Q(S) = psiT + sum lam + sum u + sum wc([I in S] - [J in S]) >= 0 on independent sets,
Phi(t) = psi0 Om(t) - sum lam - sum u_ij Om(r_ij t) >= m for t > 0, c = Phi(0) - m > 0; then m delta + c delta^2 <= psiT.
usage: python3 m1cert3.py CONF.json OUT.json [--target=0.1515] [--step=0.02] [--safety=2e-5] [--tailmargin=0.002] [--reuse]
(--reuse: take the first program's solution from OUT.lp.npz; the certificate is then rebuilt byte for byte)"""
import sys, json, time, os
from fractions import Fraction as F
import numpy as np, highspy, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conf3, ielp3
om = ielp3.om
opts = dict(a.split('=', 1) for a in sys.argv[3:] if a.startswith('--') and '=' in a)
target = float(opts.get('--target', 0.1515)); safety = float(opts.get('--safety', 2e-5))
tailmargin = float(opts.get('--tailmargin', 0.002))
pts = conf3.load(sys.argv[1]); out = sys.argv[2]; n = len(pts)
X = conf3.floats(pts)
D2, edges = conf3.graph(pts); atoms = conf3.independent_sets(n, edges); pairs, ncl = conf3.congruences(n, D2, atoms)
adj = np.zeros((n, n), dtype=bool)
for i, j in edges: adj[i, j] = adj[j, i] = True
M = np.array(atoms, dtype=np.int64)
rel = [(sum(1 << a for a in I), sum(1 << b for b in J)) for I, J in pairs]
print(f"{n} points, {len(edges)} unit pairs, {len(atoms)} independent sets, {ncl} classes, {len(rel)} relations", flush=True)
T0 = time.time()
cuts = np.arange(0.10, 0.20, 0.0025)
step = float(opts.get('--step', 0.02))
tgrid = np.concatenate([np.arange(0, 60, step), np.arange(60, 100, 0.1), np.arange(100, 600.01, 0.25)])
cache = out + '.lp.npz'
npair = [(i, j) for i in range(n) for j in range(i + 1, n) if not adj[i, j]]
if '--reuse' in sys.argv and os.path.exists(cache):
    z = np.load(cache); v, y = float(z['value']), z['y']
else:
    v, dt, info = ielp3.solve(X, tgrid=tgrid, rel=rel, M=M, adj=adj, want_dual=True, cuts=cuts)
    print(f"program value {v:.8f} ({info['status']}, {dt:.0f}s)", flush=True)
    y = info['y']; np.savez(cache, value=v, y=y)
# coarse duals: psi0 = -y0, psiT = -y1, lam, u, w; cut multipliers mu_a >= 0 give Phi(t) >= 1 - 2 sum a mu_a =: m for
# t > 0 and Phi(0) - m = sum mu_a =: c. Normalise c = 1 (divide by c).
nc = len(cuts); ycut = y[-nc:]; cmu = float(np.abs(ycut).sum())
print(f"coarse cut multipliers: sum {cmu:.6f}", flush=True)
wc = -y[2 + n + len(npair):2 + n + len(npair) + len(rel)] / cmu
G = 2 ** 40
Wi = [int(round(float(x) * G)) for x in wc]
WSa = np.zeros(len(M), dtype=np.int64)
for (I, J), w in zip(rel, Wi):
    if w: WSa[(M & I) == I] += w; WSa[(M & J) == J] -= w
Wf = WSa.astype(np.float64) / G
pr = np.sqrt(np.array([float(D2[i][j]) for i, j in npair])); m = len(npair); nv = 3 + n + m
# variables x = (psi0, psiT, mm, lam_1..n, u_1..m); Phi(t) - mm >= safety (t > 0), Phi(0) - mm = 1, -sum lam - mm >= tailmargin
# bound: delta^2 + mm delta - psiT <= 0; objective psiT - target mm (tangent at delta = target)
L0 = 3
bits = ((M[:, None] >> np.arange(n)[None, :]) & 1).astype(np.float64)
PB = bits[:, [i for i, j in npair]] * bits[:, [j for i, j in npair]]
h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', 1)
cost = np.zeros(nv); cost[1] = 1.0; cost[2] = -target
h.addVars(nv, np.full(nv, -highspy.kHighsInf), np.full(nv, highspy.kHighsInf))
h.changeColsCost(nv, np.arange(nv, dtype=np.int32), cost)
def addrows(A, lower, upper=None):
    S = sp.csr_matrix(np.asarray(A)); k = S.shape[0]
    h.addRows(k, np.asarray(lower, dtype=float), np.full(k, highspy.kHighsInf) if upper is None else np.asarray(upper, float),
              S.nnz, S.indptr[:-1].astype(np.int32), S.indices.astype(np.int32), S.data)
def phirows(t):
    t = np.asarray(t, dtype=float)
    A = np.zeros((len(t), nv)); A[:, 0] = om(t); A[:, 2] = -1.0; A[:, L0:L0 + n] = -1.0; A[:, L0 + n:] = -om(t[:, None] * pr[None, :])
    addrows(A, np.full(len(t), safety))
phirows(np.arange(0.05, 100, 0.1))
A = np.zeros((1, nv)); A[0, L0:L0 + n] = -1.0; A[0, 2] = -1.0; addrows(A, [tailmargin])
A = np.zeros((1, nv)); A[0, 0] = 1; A[0, 2] = -1; A[0, L0:] = -1; addrows(A, [1.0], [1.0])      # Phi(0) - mm = 1
for c0_ in range(0, len(M), 5000):
    k = min(5000, len(M) - c0_)
    A = np.zeros((k, nv)); A[:, 1] = 1.0; A[:, L0:L0 + n] = bits[c0_:c0_ + k]; A[:, L0 + n:] = PB[c0_:c0_ + k]
    addrows(A, -Wf[c0_:c0_ + k])
tchk = np.concatenate([np.arange(1e-3, 50, 0.005), np.arange(50, 1000, 0.02)])
Ochk = om(tchk[:, None] * pr[None, :]); O0 = om(tchk)
h.setOptionValue('solver', 'simplex')
def bound(psiT, mm): return (-mm + (mm * mm + 4 * psiT) ** 0.5) / 2
for it in range(50):
    h.run(); x = np.array(h.getSolution().col_value)
    ph = x[0] * O0 - x[L0:L0 + n].sum() - Ochk @ x[L0 + n:] - x[2]
    bad = np.flatnonzero((ph[1:-1] < safety * 0.999) & (ph[1:-1] <= ph[:-2]) & (ph[1:-1] <= ph[2:])) + 1
    print(f"polishing round {it}: {h.modelStatusToString(h.getModelStatus())}, psiT {x[1]:.8f}, m {x[2]:.6f}, bound "
          f"{bound(x[1], x[2]):.8f}, min Phi - m {ph.min():.3e}, {len(bad)} local minima below", flush=True)
    if len(bad) == 0 or ph.min() >= safety * 0.75: break
    phirows(tchk[np.unique(np.clip(np.concatenate([bad - 1, bad, bad + 1]), 0, len(tchk) - 1))])
psi0, mm, lam, u = x[0], x[2], x[L0:L0 + n], x[L0 + n:]
q = lambda z: F(int(round(float(z) * G)), G)
qd = lambda z: F(int(np.floor(float(z) * G)), G)                # rounded down onto the grid 2^-40
P0 = q(psi0); L = [q(z) for z in lam]; U = [q(z) for z in u]; Wfr = [F(w, G) for w in Wi]
def phi(t, p0, L, U):
    p = p0 * om(t) - L.sum()
    for r, w in zip(pr, U):
        if w != 0: p = p - w * om(r * t)
    return p
tv = np.concatenate([np.linspace(1e-6, 50, 2000001), np.linspace(50, 5000, 500001), np.linspace(5000, 10000, 250001)])
fl = lambda Z: np.array([float(z) for z in Z])
ph = phi(tv, float(P0), fl(L), fl(U)); mphi = ph.min()
MM = qd(mphi - safety / 2)                                   # the certified level m, a little below min Phi
print(f"min Phi on the grid {mphi:.9f} at t = {tv[ph.argmin()]:.4f}; m = {float(MM):.9f} (program m {mm:.9f})", flush=True)
Lq = [int(z * G) for z in L]; Uq = [int(z * G) for z in U]; Wq = [int(w * G) for w in Wfr]
Qa = np.zeros(len(M), dtype=np.int64)
for i in range(n):
    if Lq[i]: Qa[(M >> i) & 1 == 1] += Lq[i]
for (i, j), w in zip(npair, Uq):
    if w: Qa[((M >> i) & 1 == 1) & ((M >> j) & 1 == 1)] += w
for (I, J), w in zip(rel, Wq):
    if w: Qa[(M & I) == I] += w; Qa[(M & J) == J] -= w
qmin = int(Qa[1:].min())
PT = F(max(0, -qmin), G) + F(1, 10 ** 9); PT = F(-((-PT * G) // 1), G)
C0 = P0 - sum(L) - sum(U) - MM
if C0 <= 0: raise SystemExit('Phi(0) - m <= 0')
B = (-float(MM) + (float(MM) ** 2 + 4 * float(C0) * float(PT)) ** 0.5) / (2 * float(C0))
print(f"psiT = {float(PT):.10f}, m = {float(MM):.8f}, c = Phi(0) - m = {float(C0):.8f}: bound {B:.10f} = 1/{1 / B:.5f} "
      f"(program {v:.8f})", flush=True)
cert = dict(side='M1R3', scale='1/sqrt2', points=[[str(c) for c in p] for p in pts], psi0=str(P0), psiT=str(PT), m=str(MM),
            lam=[str(z) for z in L], u=[[i + 1, j + 1, str(w)] for (i, j), w in zip(npair, U) if w],
            wc=[[[a + 1 for a in I], [b + 1 for b in J], str(w)] for (I, J), w in zip(pairs, Wfr) if w])
json.dump(cert, open(out, 'w'))
print(f"wrote {out} ({time.time() - T0:.0f}s)")
