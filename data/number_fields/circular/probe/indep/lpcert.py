"""Exact optimum of a 3-variable LP  min c.z  s.t.  A z <= b  (Fractions), found with a floating-point solver and then
certified exactly: a vertex given by 3 tight constraints that is primal feasible (exact check against all rows) and
has nonnegative exact KKT multipliers (c + sum m_i a_i = 0)."""
from fractions import Fraction as Fr
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
from lp3x import solve3


def exact_lp(c, rows):
    A = np.array([[float(r[0]), float(r[1]), float(r[2])] for r in rows])
    b = np.array([float(r[3]) for r in rows])
    res = linprog([float(t) for t in c], A_ub=A, b_ub=b, bounds=[(None, None)] * 3, method="highs")
    if res.status == 2:
        return None  # infeasible (float); caller may double-check
    assert res.status == 0, res.message
    z = res.x
    slack = b - A @ z
    order = np.argsort(slack)
    cand = [int(i) for i in order[:12]]
    for tri in combinations(cand, 3):
        v = solve3([rows[i] for i in tri])
        if v is None:
            continue
        if not all(r[0] * v[0] + r[1] * v[1] + r[2] * v[2] <= r[3] for r in rows):
            continue
        # multipliers: sum m_i a_i = -c
        M = [[rows[i][t] for i in tri] for t in range(3)]
        rhs = [-Fr(t) for t in c]
        mrow = [M[t] + [rhs[t]] for t in range(3)]
        m = solve3(mrow)
        if m is None or any(x < 0 for x in m):
            continue
        val = sum(Fr(c[t]) * v[t] for t in range(3))
        return val, v, tri, m
    raise RuntimeError("could not certify")
