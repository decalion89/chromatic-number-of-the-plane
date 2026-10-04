"""Exact 3-variable LP:   maximise t  subject to  t <= a_i x + b_i y + c_i   (i = 1..n).

A floating-point solve (scipy HiGHS) proposes an optimal vertex; the result is returned only with an exact
certificate: an exactly feasible point (x*, y*, t*) and exact multipliers lam_i >= 0 on three constraints with
sum lam_i (-a_i, -b_i, 1) = (0, 0, 1); then for every feasible (x, y, t),
    t = sum lam_i t <= sum lam_i (a_i x + b_i y + c_i) = sum lam_i c_i = t*      (weak duality),
so t* is the exact optimum.  If no certificate is found among the near-active constraints, all triples of
constraints are tried exactly (slow but complete)."""
from fractions import Fraction as Fr
from itertools import combinations


def solve3(M, v):
    (a, b, c), (d, e, f), (g, h, i) = M
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if det == 0:
        return None
    x = (v[0] * (e * i - f * h) - b * (v[1] * i - f * v[2]) + c * (v[1] * h - e * v[2])) / det
    y = (a * (v[1] * i - f * v[2]) - v[0] * (d * i - f * g) + c * (d * v[2] - v[1] * g)) / det
    z = (a * (e * v[2] - v[1] * h) - b * (d * v[2] - v[1] * g) + v[0] * (d * h - e * g)) / det
    return (x, y, z)


def certify(cons, tri):
    """cons: list of (a, b, c). tri: 3 indices. Returns (t, x, y, lam) if the vertex is feasible and dual-optimal."""
    rows = [(-cons[i][0], -cons[i][1], Fr(1)) for i in tri]
    rhs = [cons[i][2] for i in tri]
    sol = solve3(rows, rhs)
    if sol is None:
        return None
    x, y, t = sol
    for (a, b, c) in cons:
        if t > a * x + b * y + c:
            return None
    MT = [[rows[0][j], rows[1][j], rows[2][j]] for j in range(3)]
    lam = solve3(MT, (Fr(0), Fr(0), Fr(1)))
    if lam is None or any(l < 0 for l in lam):
        return None
    return (t, x, y, dict(zip(tri, lam)))


def lp_max(cons):
    import numpy as np
    from scipy.optimize import linprog
    A = np.array([[-float(a), -float(b), 1.0] for (a, b, c) in cons])
    bv = np.array([float(c) for (a, b, c) in cons])
    res = linprog(c=[0, 0, -1], A_ub=A, b_ub=bv, bounds=[(None, None)] * 3, method="highs")
    if res.status == 0:
        slack = bv - A @ res.x
        order = [int(i) for i in np.argsort(slack)[:12]]
        for tri in combinations(order, 3):
            c = certify(cons, tri)
            if c is not None:
                return c
    # complete fallback
    best = None
    for tri in combinations(range(len(cons)), 3):
        c = certify(cons, tri)
        if c is not None and (best is None or c[0] > best[0]):
            best = c
    return best
