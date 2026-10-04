"""Exact LP in three variables by vertex enumeration (the feasible sets used here are pointed polyhedra; the
objective is bounded below on them).  Constraints are (a1, a2, a3, b) meaning a1 z1 + a2 z2 + a3 z3 <= b."""
from fractions import Fraction as Fr
from itertools import combinations


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def solve3(rows):
    A = [r[:3] for r in rows]
    b = [r[3] for r in rows]
    D = det3(A)
    if D == 0:
        return None
    out = []
    for c in range(3):
        M = [list(A[i]) for i in range(3)]
        for i in range(3):
            M[i][c] = b[i]
        out.append(det3(M) / D)
    return tuple(out)


def lp_min(obj, cons):
    """min obj.z over {z : cons}; returns (value, z, tight-index-triple) or None if infeasible"""
    best = None
    seen = set()
    for tri in combinations(range(len(cons)), 3):
        z = solve3([cons[i] for i in tri])
        if z is None or z in seen:
            continue
        seen.add(z)
        if all(c[0] * z[0] + c[1] * z[1] + c[2] * z[2] <= c[3] for c in cons):
            val = obj[0] * z[0] + obj[1] * z[1] + obj[2] * z[2]
            if best is None or val < best[0]:
                best = (val, z, tri)
    return best


def farkas_check(cons, mult, idx, target):
    """check sum mult_i * cons[idx_i] reads 0*z1 + 0*z2 - z3 <= -target, i.e. z3 >= target"""
    tot = [Fr(0)] * 4
    for m, i in zip(mult, idx):
        assert m >= 0
        for t in range(4):
            tot[t] += m * cons[i][t]
    return tot[0] == 0 and tot[1] == 0 and tot[2] == -1 and tot[3] == -target
