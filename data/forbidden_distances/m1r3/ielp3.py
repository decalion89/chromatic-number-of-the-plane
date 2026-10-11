"""Inclusion-exclusion program of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki in R^3 (floating point; exploration).
X: (n, 3) array. Radial positive-definite kernel in R^3: Omega_3(s) = sin(s)/s. Congruence = equal squared-distance
matrices under a bijection (any distance-preserving bijection of finite sets extends to an isometry of R^3)."""
import time
from collections import defaultdict
import numpy as np, scipy.sparse as sp
import highspy

TOL = 1e-9
def om(s): return np.sinc(np.asarray(s, dtype=float) / np.pi)

def unit_graph(X):
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2); return D, np.abs(D - 1) < TOL

def indep_sets(n, adj):
    nb = [sum(1 << u for u in range(n) if adj[v, u]) for v in range(n)]
    out = []
    def rec(start, cur, mask):
        out.append(cur)
        for v in range(start, n):
            if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nb[v])
    rec(0, 0, 0)
    return np.array(out, dtype=np.int64)

def congruences(D, M, nd=7, minsize=2):
    """relations (mask_I, mask_J) between congruent independent sets with >= minsize points"""
    n = len(D); Dr = np.round(D, nd)
    def sig(S): return tuple(sorted(Dr[a, b] for a in S for b in S if a < b))
    def iso(S, T):
        k = len(S); used = [False] * k; img = [None] * k
        def bt(p):
            if p == k: return True
            for c in range(k):
                if used[c]: continue
                if all(Dr[S[q], S[p]] == Dr[T[img[q]], T[c]] for q in range(p)):
                    used[c] = True; img[p] = c
                    if bt(p + 1): return True
                    used[c] = False
            return False
        return bt(0)
    groups = defaultdict(list)
    for m in M:
        m = int(m)
        S = tuple(i for i in range(n) if m >> i & 1)
        if len(S) < minsize: continue
        groups[(len(S), sig(S))].append((m, S))
    rel = []; ncl = 0
    for mem in groups.values():
        reps = []
        for m, S in mem:
            for rm, rS in reps:
                if iso(rS, S): rel.append((rm, m)); break
            else: reps.append((m, S))
        ncl += len(reps)
    return rel, ncl

def solve(X, tgrid=None, rel=None, M=None, adj=None, solver='ipm', threads=1, crossover='off', want_dual=False, verbose=False,
          cuts=()):
    """cuts: values a for the tangent cuts kappa(0) >= 2 a delta - a^2 of kappa(0) >= delta^2 (fhat(0) = delta^2)"""
    n = len(X); D, adj0 = unit_graph(X)
    if adj is None: adj = adj0
    if M is None: M = indep_sets(n, adj)
    if rel is None: rel, _ = congruences(D, M, minsize=3)
    t = np.concatenate([np.arange(0, 100, 0.1), np.arange(100, 400.01, 0.25)]) if tgrid is None else np.asarray(tgrid)
    nk = len(t); na = len(M)
    R, C, V = [], [], []; rhs = []
    def add(kcoef=None, acols=None, avals=None, b=0.0):
        r = len(rhs)
        if kcoef is not None:
            nz = np.flatnonzero(kcoef); R.append(np.full(len(nz), r)); C.append(nz); V.append(kcoef[nz])
        if acols is not None:
            acols = np.asarray(acols); R.append(np.full(len(acols), r)); C.append(nk + acols)
            V.append(np.broadcast_to(np.asarray(avals, dtype=float), len(acols)).astype(float))
        rhs.append(b)
    add(kcoef=om(t))
    add(acols=np.arange(na), avals=1.0, b=1.0)
    for i in range(n): add(kcoef=-np.ones(nk), acols=np.flatnonzero((M >> i) & 1), avals=1.0)
    for i in range(n):
        for j in range(i + 1, n):
            if not adj[i, j]: add(kcoef=-om(t * D[i, j]), acols=np.flatnonzero(((M >> i) & 1) & ((M >> j) & 1)), avals=1.0)
    for I, J in rel:
        a1 = np.flatnonzero((M & I) == I); a2 = np.flatnonzero((M & J) == J)
        add(acols=np.concatenate([a1, a2]), avals=np.concatenate([np.ones(len(a1)), -np.ones(len(a2))]))
    neq = len(rhs)
    assert t[0] == 0
    for a in cuts:
        kc = np.full(nk, -2.0 * a); kc[0] += 1.0; add(kcoef=kc, b=-a * a)
    A = sp.csc_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(len(rhs), nk + na))
    up = np.array(rhs); up[neq:] = highspy.kHighsInf
    cost = np.concatenate([-np.ones(nk), np.zeros(na)])
    h = highspy.Highs(); h.setOptionValue('output_flag', verbose)
    h.setOptionValue('solver', solver); h.setOptionValue('threads', threads); h.setOptionValue('run_crossover', crossover)
    lp = highspy.HighsLp(); lp.num_col_ = A.shape[1]; lp.num_row_ = A.shape[0]; lp.col_cost_ = cost
    lp.col_lower_ = np.zeros(A.shape[1]); lp.col_upper_ = np.full(A.shape[1], highspy.kHighsInf)
    lp.row_lower_ = np.array(rhs); lp.row_upper_ = up
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    h.passModel(lp); t0 = time.time(); h.run(); dt = time.time() - t0
    v = -h.getInfo().objective_function_value
    info = dict(status=h.modelStatusToString(h.getModelStatus()), atoms=na, rel=len(rel), edges=int(adj.sum()) // 2)
    if want_dual:
        sol = h.getSolution(); info.update(y=np.array(sol.row_dual), x=np.array(sol.col_value), t=t, nk=nk)
    return v, dt, info
