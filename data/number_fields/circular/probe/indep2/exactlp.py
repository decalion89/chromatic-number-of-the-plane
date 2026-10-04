# Exact 3-variable LP: minimise the 3rd variable s subject to A z <= b (z=(z1,z2,s)).
# Two methods: brute-force vertex enumeration (small instances), and float LP (scipy HiGHS)
# followed by an exact primal/dual certificate on a basis of 3 tight constraints.
from fractions import Fraction as Fr
import itertools

def solve3(M, r):
    """Solve 3x3 system M z = r exactly; return None if singular."""
    a = [list(map(Fr, M[i])) + [Fr(r[i])] for i in range(3)]
    for col in range(3):
        piv = next((i for i in range(col, 3) if a[i][col] != 0), None)
        if piv is None:
            return None
        a[col], a[piv] = a[piv], a[col]
        for i in range(3):
            if i != col and a[i][col] != 0:
                f = a[i][col] / a[col][col]
                a[i] = [a[i][t] - f * a[col][t] for t in range(4)]
    return [a[i][3] / a[i][i] for i in range(3)]

def brute_min_s(A, b):
    best = None; arg = None
    for T in itertools.combinations(range(len(A)), 3):
        z = solve3([A[t] for t in T], [b[t] for t in T])
        if z is None: continue
        if best is not None and z[2] >= best: continue
        if all(sum(A[i][t] * z[t] for t in range(3)) <= b[i] for i in range(len(A))):
            best = z[2]; arg = z
    return best, arg

def certified_min_s(A, b):
    """Return (min s, z, dual) with an exact certificate, or (None, farkas) if infeasible."""
    import numpy as np
    from scipy.optimize import linprog
    Af = np.array([[float(x) for x in row] for row in A]); bf = np.array([float(x) for x in b])
    res = linprog(c=[0, 0, 1], A_ub=Af, b_ub=bf, bounds=[(None, None)] * 3, method="highs")
    if res.status == 2:
        return None, None, None
    assert res.status == 0, res.message
    z0 = res.x
    slack = bf - Af @ z0
    order = sorted(range(len(A)), key=lambda i: abs(slack[i]))
    cand = order[:12]
    for T in itertools.combinations(cand, 3):
        M = [A[t] for t in T]
        z = solve3(M, [b[t] for t in T])
        if z is None: continue
        if not all(sum(A[i][t] * z[t] for t in range(3)) <= b[i] for i in range(len(A))):
            continue
        # dual: find y >= 0 with sum_t y_t A[t] = -(0,0,1)  (KKT for min s with A z <= b).
        # Then for every feasible z: -s = sum y_t A[t].z <= sum y_t b_t, i.e. s >= -sum y_t b_t.
        MT = [[M[r][c] for r in range(3)] for c in range(3)]
        y = solve3(MT, [0, 0, -1])
        if y is None or any(v < 0 for v in y):
            continue
        val = -sum(y[i] * b[T[i]] for i in range(3))
        assert val == z[2]
        return z[2], z, dict(zip(T, y))
    raise RuntimeError("no exact certificate found among near-tight constraints")
