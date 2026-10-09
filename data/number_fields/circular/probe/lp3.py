"""Exact 3-variable linear programming with certificates.

lp3_max(obj, rows): maximise obj . z over z in Q^3 subject to a . z <= b for (a1, a2, a3, b) in rows.
Returns (value, z*, active_rows, multipliers) where z* is an exactly feasible vertex, active_rows three
linearly independent rows tight at z*, and multipliers >= 0 with sum m_i a_i = obj; by weak duality
obj . z <= sum m_i b_i = obj . z* for every feasible z, so the value is exact and certified.
The search for the vertex uses a float LP (scipy/HiGHS) and then exact rational arithmetic; if that fails,
exact enumeration of all vertices is used (the LP must be bounded with a pointed feasible region).
check_certificate() re-verifies a returned certificate from scratch.
"""
from fractions import Fraction as F
from itertools import combinations


def det3(M):
    (a, b, c), (d, e, f), (g, h, i) = M
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def solve3(M, v):
    D = det3(M)
    if D == 0:
        return None
    out = []
    for col in range(3):
        Mc = [list(row) for row in M]
        for rr in range(3):
            Mc[rr][col] = v[rr]
        out.append(det3(Mc) / D)
    return tuple(out)


def _cert(obj, rows, tri):
    M = [rows[i][:3] for i in tri]
    v = [rows[i][3] for i in tri]
    z = solve3(M, v)
    if z is None:
        return None
    if not all(r[0] * z[0] + r[1] * z[1] + r[2] * z[2] <= r[3] for r in rows):
        return None
    MT = [[M[0][c], M[1][c], M[2][c]] for c in range(3)]
    m = solve3(MT, tuple(obj))
    if m is None or any(x < 0 for x in m):
        return None
    val = obj[0] * z[0] + obj[1] * z[1] + obj[2] * z[2]
    return (val, z, tri, m)


def lp3_max(obj, rows, brute_ok=True):
    import numpy as np
    from scipy.optimize import linprog
    obj = tuple(F(x) for x in obj)
    rows = [tuple(F(x) for x in r) for r in rows]
    A = np.array([[float(r[0]), float(r[1]), float(r[2])] for r in rows])
    b = np.array([float(r[3]) for r in rows])
    res = linprog(c=[-float(x) for x in obj], A_ub=A, b_ub=b, bounds=[(None, None)] * 3, method="highs")
    if res.status == 2:
        return None          # infeasible (float); verified separately if needed
    if res.status == 0:
        slack = b - A @ res.x
        order = [i for i in np.argsort(slack)][:10]
        for tri in combinations(order, 3):
            c = _cert(obj, rows, tri)
            if c is not None:
                return c
    if not brute_ok:
        raise RuntimeError("no certificate found")
    best = None
    for tri in combinations(range(len(rows)), 3):
        M = [rows[i][:3] for i in tri]
        v = [rows[i][3] for i in tri]
        z = solve3(M, v)
        if z is None:
            continue
        if all(r[0] * z[0] + r[1] * z[1] + r[2] * z[2] <= r[3] for r in rows):
            val = obj[0] * z[0] + obj[1] * z[1] + obj[2] * z[2]
            if best is None or val > best[0]:
                best = (val, z, tri)
    if best is None:
        return None
    # find a dual certificate among all tight triples at the best vertex
    z = best[1]
    tight = [i for i, r in enumerate(rows) if r[0] * z[0] + r[1] * z[1] + r[2] * z[2] == r[3]]
    for tri in combinations(tight, 3):
        c = _cert(obj, rows, tri)
        if c is not None:
            return c
    raise RuntimeError("optimum without certificate")


def check_certificate(obj, rows, cert):
    val, z, tri, m = cert
    obj = tuple(F(x) for x in obj)
    rows = [tuple(F(x) for x in r) for r in rows]
    assert all(r[0] * z[0] + r[1] * z[1] + r[2] * z[2] <= r[3] for r in rows), "primal infeasible"
    assert all(x >= 0 for x in m), "negative multiplier"
    for c in range(3):
        assert sum(m[t] * rows[tri[t]][c] for t in range(3)) == obj[c], "dual equation"
    assert sum(m[t] * rows[tri[t]][3] for t in range(3)) == val, "duality gap"
    assert obj[0] * z[0] + obj[1] * z[1] + obj[2] * z[2] == val
    return True
