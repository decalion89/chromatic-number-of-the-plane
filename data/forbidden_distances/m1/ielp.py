"""One-distance inclusion-exclusion program of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki for an arbitrary finite
point set X (floating point; exploration only, nothing here is a proof).
maximize sum kappa  s.t.  kappa, p >= 0;  sum kappa(t) J0(t) = 0;  sum_S p(S) = 1;  P(i in S) = sum kappa;
P(i, j in S) = sum kappa(t) J0(r_ij t) for non-unit pairs;  P(I in S) = P(J in S) for congruent independent I, J.
S runs over the independent sets of the unit-distance graph of X."""
import time
import numpy as np, scipy.sparse as sp
from scipy.special import j0
import highspy

TOL = 1e-9

def unit_graph(X):
    D = np.abs(X[:, None] - X[None, :]); return D, np.abs(D - 1) < TOL

def indep_sets(n, adj):
    nb = [sum(1 << u for u in range(n) if adj[v, u]) for v in range(n)]
    out = []
    def rec(start, cur, mask):
        out.append(cur)
        for v in range(start, n):
            if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nb[v])
    rec(0, 0, 0)
    return np.array(out, dtype=np.int64)

def canon(P, nd=6):
    """isometry-invariant canonical form of the point set P (complex array, len >= 2)"""
    best = None
    k = len(P)
    for a in range(k):
        for b in range(k):
            if a == b: continue
            v = P[b] - P[a]; w = (P - P[a]) * np.conj(v) / abs(v)
            for z in (w, np.conj(w)):
                key = tuple(sorted((round(x.real, nd) + 0.0, round(x.imag, nd) + 0.0) for x in z))
                if best is None or key < best: best = key
    return best

def congruences(X, M):
    """relations (mask_I, mask_J): one per non-representative member of each congruence class of independent sets
    with at least two points"""
    classes = {}
    for m in M:
        m = int(m)
        if m & (m - 1) == 0: continue                        # empty set and singletons
        idx = [i for i in range(len(X)) if m >> i & 1]
        classes.setdefault(canon(X[idx]), []).append(m)
    rel = []
    for ms in classes.values():
        for m in ms[1:]: rel.append((ms[0], m))
    return rel, len(classes)

def build(X, tmax=100.0, step=0.1, rel=None, M=None, adj=None, tgrid=None):
    """constraint matrix (scipy csc), right-hand side, t grid, number of atoms, relations (tgrid, if given, replaces
    the grid 0, step, ..., tmax)"""
    n = len(X); D, adj0 = unit_graph(X)
    if adj is None: adj = adj0
    if M is None: M = indep_sets(n, adj)
    if rel is None: rel, _ = congruences(X, M)
    M = np.asarray(M, dtype=np.int64)
    t = np.arange(0, tmax + 1e-9, step) if tgrid is None else np.asarray(tgrid, dtype=float); nk = len(t); na = len(M)
    R, C, V = [], [], []                                  # lists of numpy arrays (rows, columns, values)
    rhs = []
    def add(kcoef=None, acols=None, avals=None, b=0.0):
        r = len(rhs)
        if kcoef is not None:
            nz = np.flatnonzero(kcoef)
            R.append(np.full(len(nz), r, dtype=np.int32)); C.append(nz.astype(np.int32)); V.append(kcoef[nz])
        if acols is not None:
            acols = np.asarray(acols)
            R.append(np.full(len(acols), r, dtype=np.int32)); C.append((nk + acols).astype(np.int32))
            V.append(np.broadcast_to(np.asarray(avals, dtype=float), len(acols)).astype(float))
        rhs.append(b)
    add(kcoef=j0(t))
    add(acols=np.arange(na), avals=1.0, b=1.0)
    for i in range(n): add(kcoef=-np.ones(nk), acols=np.flatnonzero((M >> i) & 1), avals=1.0)
    for i in range(n):
        for j in range(i + 1, n):
            if not adj[i, j]:
                add(kcoef=-j0(t * D[i, j]), acols=np.flatnonzero(((M >> i) & 1) & ((M >> j) & 1)), avals=1.0)
    for I, J in rel:
        a1 = np.flatnonzero((M & I) == I); a2 = np.flatnonzero((M & J) == J)
        add(acols=np.concatenate([a1, a2]), avals=np.concatenate([np.ones(len(a1)), -np.ones(len(a2))]))
    R = np.concatenate(R); C = np.concatenate(C); V = np.concatenate(V)
    A = sp.csc_matrix((V, (R, C)), shape=(len(rhs), nk + na))
    del R, C, V
    return A, np.array(rhs), t, nk, na, rel

def solve(X, tmax=100.0, step=0.1, solver='ipm', threads=1, verbose=False, rel=None, M=None, want_dual=False, adj=None,
          acost=None, crossover='on', tgrid=None):
    """value of the program; solver 'ipm' or 'simplex' (HiGHS) or 'clarabel' (threads = its thread count). With
    want_dual, info has y (row duals in the HiGHS sign convention: psi0 = -y[0], psiT = -y[1], ...) and x"""
    A, rhs, t, nk, na, rel = build(X, tmax, step, rel, M, adj, tgrid)
    cost = np.concatenate([-np.ones(nk), np.zeros(na) if acost is None else np.asarray(acost, dtype=float)])
    if solver == 'clarabel':
        import clarabel
        m, nv = A.shape
        Ab = sp.vstack([A, -sp.identity(nv, format='csc')]).tocsc()
        b = np.concatenate([rhs, np.zeros(nv)])
        st = clarabel.DefaultSettings(); st.verbose = verbose
        st.tol_gap_abs = 1e-9; st.tol_gap_rel = 1e-9; st.tol_feas = 1e-9
        st.direct_solve_method = 'faer'; st.max_threads = threads
        t0 = time.time()
        s = clarabel.DefaultSolver(sp.csc_matrix((nv, nv)), cost, Ab, b, [clarabel.ZeroConeT(m), clarabel.NonnegativeConeT(nv)], st)
        sol = s.solve(); dt = time.time() - t0
        v = -sol.obj_val
        info = dict(status=str(sol.status), rows=m, atoms=na, rel=len(rel))
        if want_dual:
            z = np.array(sol.z)
            info.update(y=-z[:m], x=np.array(sol.x), t=t, nk=nk)
        return v, dt, info
    h = highspy.Highs(); h.setOptionValue('output_flag', verbose)
    h.setOptionValue('solver', solver); h.setOptionValue('threads', threads); h.setOptionValue('run_crossover', crossover)
    lp = highspy.HighsLp()
    lp.num_col_ = A.shape[1]; lp.num_row_ = A.shape[0]
    lp.col_cost_ = cost
    lp.col_lower_ = np.zeros(A.shape[1]); lp.col_upper_ = np.full(A.shape[1], highspy.kHighsInf)
    lp.row_lower_ = rhs; lp.row_upper_ = rhs
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data
    h.passModel(lp)
    t0 = time.time(); h.run(); dt = time.time() - t0
    st = h.getModelStatus()
    v = -h.getInfo().objective_function_value
    info = dict(status=h.modelStatusToString(st), rows=A.shape[0], atoms=na, rel=len(rel))
    if want_dual:
        sol = h.getSolution()
        info.update(y=np.array(sol.row_dual), x=np.array(sol.col_value), t=t, nk=nk)
    return v, dt, info
