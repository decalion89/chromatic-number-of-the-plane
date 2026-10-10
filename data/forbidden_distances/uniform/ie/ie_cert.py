"""Rational certificates from the floating-point duals of ie_lp.py (check_ie.py then proves them).

A-side: round the multipliers to multiples of 2^-40, raise psiT so that every atom inequality holds (and psiT >= zeta
with the hit-set row), and scale everything so that Phi (minus the hit-set term) is at least 1 + safety on a fine grid.
B-side (--B): put the multipliers on the grid 2^-36, choose psiT and the unit-pair penalty pi by an exact integer check
over all 2^23 subsets, and scale so that Phi_B >= 1 + safety on a fine grid.
usage: python3 ie_cert.py DUAL.npz OUT.json [--B] [--hit=certB.json] [--safety=1e-4]"""
import sys, os, json
from fractions import Fraction as F
import numpy as np
from scipy.special import j0
from x23 import X
opts = dict(a.split('=', 1) for a in sys.argv[3:] if a.startswith('--') and '=' in a)
Z = np.load(sys.argv[1], allow_pickle=True); out = sys.argv[2]; safety = float(opts.get('--safety', 1e-4))
BSIDE = '--B' in sys.argv
MU = 0.4028
y = Z['yeq']; z = Z['zub']; pairs = [tuple(int(a) for a in p) for p in Z['pairs']]
cong = [tuple(c) for c in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cong23.json')))]
n = len(X); npair = len(pairs); Dm = np.abs(X[:, None] - X[None, :]); adj = np.abs(Dm - 1) < 1e-9
edges = [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i, j]]
assert len(y) == 2 + n + npair + len(cong)
psi0, psiT = -y[0], -y[1]; lam = -y[2:2 + n]; u = -y[2 + n:2 + n + npair]; wc = -y[2 + n + npair:]
masks = [(sum(1 << (a - 1) for a in I), sum(1 << (b - 1) for b in J)) for I, J in cong]
def phi_grid(t, extra=0.0):
    p = psi0 * j0(t) - lam.sum() + extra
    for (i, j), w in zip(pairs, u):
        if abs(w) > 1e-14: p = p - w * j0(t * Dm[i, j])
    return p
if BSIDE:
    t = np.concatenate([np.linspace(0, 50, 500001), np.linspace(50, 3000, 300001)])
    m = phi_grid(t).min(); s = 1 / (m - safety)
    print(f"min Phi_B on the grid {m:.8f}; scale {s:.8f}")
    G = 2 ** 36; qz = lambda v: int(round(float(v) * s * G))
    P0, Lm = qz(psi0), [qz(v) for v in lam]
    Uu = [(i, j, qz(w)) for (i, j), w in zip(pairs, u) if qz(w)]
    Wc = [(k, qz(w)) for k, w in enumerate(wc) if qz(w)]
    T = np.arange(1 << n, dtype=np.int64); Q = np.zeros(1 << n, dtype=np.int64)
    for i in range(n): Q += Lm[i] * ((T >> i) & 1)
    for i, j, w in Uu: Q += w * (((T >> i) & 1) & ((T >> j) & 1))
    for k, w in Wc:
        I, J = masks[k]; Q += w * (((T & I) == I).astype(np.int64) - ((T & J) == J).astype(np.int64))
    E = np.zeros(1 << n, dtype=np.int64)
    for i, j in edges: E += ((T >> i) & 1) & ((T >> j) & 1)
    ind = E == 0
    PT = int(-Q[ind].min()) + 1
    Q += PT
    # one penalty per unit pair: minimise their sum subject to Q(T) + sum_{e in T} pi_e >= 0 for every T (cutting
    # planes over the subsets with Q(T) < 0, in floating point), then round up and re-check exactly
    from scipy.optimize import linprog
    Em = np.zeros((len(edges), 1 << n), dtype=bool)
    for k, (i, j) in enumerate(edges): Em[k] = (((T >> i) & 1) & ((T >> j) & 1)).astype(bool)
    bad = np.flatnonzero(Q < 0); cuts = []; pi = np.zeros(len(edges))
    while True:
        slack = Q[bad] / G + Em[:, bad].T.astype(float) @ pi
        viol = np.flatnonzero(slack < -1e-12)
        if len(viol) == 0: break
        cuts = sorted(set(cuts) | set(bad[viol[np.argsort(slack[viol])[:2000]]].tolist()))
        r = linprog(np.ones(len(edges)), A_ub=-Em[:, cuts].T.astype(float), b_ub=Q[cuts] / G, bounds=(0, None),
                    method='highs')
        pi = r.x
    PI = [int(np.ceil(x * G)) + 4 for x in pi]
    Qe = Q.copy()
    for k in range(len(edges)): Qe += PI[k] * Em[k]
    assert int(Qe.min()) >= 0, 'rounding the penalties left a negative subset'
    del Qe, Em
    LBn = P0 + sum(PI)
    print(f"psiT_B = {PT / G:.8f}, per-pair penalties sum {sum(PI) / G:.4f}, L_B = psi0 + sum pi = {LBn / G:.4f}")
    cert = dict(side='B', psi0=str(F(P0, G)), psiT=str(F(PT, G)),
                pi=[[i + 1, j + 1, str(F(v, G))] for (i, j), v in zip(edges, PI)], lam=[str(F(v, G)) for v in Lm],
                u=[[i + 1, j + 1, str(F(w, G))] for i, j, w in Uu],
                wc=[[list(cong[k][0]), list(cong[k][1]), str(F(w, G))] for k, w in Wc])
else:
    D = float(Z['D']); d0 = float(Z['d0'])
    q = lambda v: F(round(float(v) * 2 ** 40), 2 ** 40)
    eta, tau = max(0.0, -z[0]), max(0.0, -z[1]); zeta = max(0.0, -z[2]) if len(z) > 2 else 0.0
    psiTB = LB = F(0)
    if zeta:
        CB = json.load(open(opts['--hit']))
        psiTB = F(CB['psiT'])
        LB = max(F(0), F(CB['psi0']) + (sum(F(v) for _, _, v in CB['pi']) if isinstance(CB['pi'], list)
                                        else len(edges) * F(CB['pi'])))
        assert float(LB) <= float(Z['LB']) + 1e-9, 'the program assumed a smaller L_B'
    nb = [sum(1 << v for v in range(n) if adj[i, v]) for i in range(n)]
    atoms = []
    def rec(start, cur, mask):
        atoms.append(cur)
        for v in range(start, n):
            if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nb[v])
    rec(0, 0, 0)
    A = np.array(atoms, dtype=np.int64); Qa = np.zeros(len(A))
    for i in range(n): Qa += float(q(lam[i])) * ((A >> i) & 1)
    for (i, j), w in zip(pairs, u): Qa += float(q(w)) * (((A >> i) & 1) & ((A >> j) & 1))
    for (I, J), w in zip(masks, wc): Qa += float(q(w)) * (((A & I) == I).astype(float) - ((A & J) == J).astype(float))
    PT = q(max(psiT, zeta, -Qa[1:].min()) + 1e-9)
    t = np.concatenate([np.geomspace(1e-6 / D, 1e3 / D, 4000), np.linspace(1e3 / D, 50, 400000),
                        np.linspace(50, 4000, 400000)])
    e = np.minimum(MU, np.sqrt(2 / (np.pi * D * t)))
    extra = -eta * e + 2 * d0 * tau
    if zeta:
        H = n + 2 * sum(j0(t * Dm[i, j]) for i in range(n) for j in range(i + 1, n))
        extra = extra - zeta * float(LB) * e * (n * n - H)
    phi = phi_grid(t, extra)
    phi0 = psi0 - lam.sum() - u.sum() + eta + tau * (2 * d0 - 1)
    m = min(phi.min(), phi0); s = 1 / (m - safety) if m - safety < 1 else 1.0
    print(f"min Phi on the grid {phi.min():.8f} (t = {t[phi.argmin()]:.6g}), Phi0 {phi0:.8f}, zeta {zeta:.6f}; scale {s:.8f}")
    S = q(s); d0F = F(d0).limit_denominator(10 ** 6)
    cert = dict(side='A', D=str(F(D).limit_denominator(1)), d0=str(d0F), mu='0.4028', psi0=str(q(psi0) * S),
                psiT=str(PT * S), eta=str(q(eta) * S), tau=str(q(tau) * S), zeta=str(q(zeta) * S),
                lam=[str(q(x) * S) for x in lam], u=[[i + 1, j + 1, str(q(w) * S)] for (i, j), w in zip(pairs, u) if q(w)],
                wc=[[list(I), list(J), str(q(w) * S)] for (I, J), w in zip(cong, wc) if q(w)])
    if zeta: cert['B'] = opts['--hit']
    Bd = (PT + q(tau) * d0F ** 2 - q(zeta) * (1 - psiTB)) * S
    print(f"bound {float(Bd):.8f} = 1/{1 / float(Bd):.4f} (program {float(Z['value']):.6f})")
json.dump(cert, open(out, 'w'))
print("wrote", out)
