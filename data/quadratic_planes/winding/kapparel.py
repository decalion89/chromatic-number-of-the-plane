"""kapparel.py d D [A] [tlimit]: Lemma W test in relation space.  f_u in [A, 1-A] for u in U_D (one per +- pair) and
<r, f> integer for every integer relation r (sum r_u v_u = 0); the relation lattice is computed exactly and
LLL-reduced (augmented-lattice trick).  Feasible iff a character xi of Z U_D has xi(U_D) in [A, 1-A]."""
import sys, time
from fractions import Fraction as Fr
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
from ortools.linear_solver import pywraplp
from kappaD import units, halve

def relation_basis(V):
    n = len(V); W = 10 ** 6
    rows = [[1 if j == i else 0 for j in range(n)] + [W * c for c in V[i]] for i in range(n)]
    M = DomainMatrix([[ZZ(x) for x in r] for r in rows], (n, n + 4), ZZ).lll()
    R = [list(map(int, M.to_Matrix().row(i))) for i in range(n)]
    ker = [r[:n] for r in R if all(x == 0 for x in r[n:])]
    for r in ker:
        assert all(sum(r[i] * V[i][j] for i in range(n)) == 0 for j in range(4))
    return ker

def solve_rel(n, ker, A, tlimit):
    s = pywraplp.Solver.CreateSolver("SCIP")
    f = [s.NumVar(float(A), float(1 - A), f"f{i}") for i in range(n)]
    for j, r in enumerate(ker):
        lo = sum(min(0, c) for c in r) * float(1 - A) + sum(max(0, c) for c in r) * float(A)
        hi = sum(max(0, c) for c in r) * float(1 - A) + sum(min(0, c) for c in r) * float(A)
        import math
        z = s.IntVar(math.floor(lo) - 1, math.ceil(hi) + 1, f"z{j}")
        s.Add(sum(c * f[i] for i, c in enumerate(r) if c) == z)
    s.SetTimeLimit(int(tlimit * 1000))
    st = s.Solve()
    name = {pywraplp.Solver.OPTIMAL: "FEASIBLE", pywraplp.Solver.FEASIBLE: "FEASIBLE",
            pywraplp.Solver.INFEASIBLE: "INFEASIBLE"}.get(st, f"UNKNOWN(status {st})")
    return name, ([x.solution_value() for x in f] if name == "FEASIBLE" else None)

if __name__ == "__main__":
    d, D = int(sys.argv[1]), int(sys.argv[2])
    A = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1, 3)
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    t0 = time.time()
    V = halve(units(d, D)); n = len(V)
    ker = relation_basis(V)
    mx = max((sum(abs(c) for c in r) for r in ker), default=0)
    name, f = solve_rel(n, ker, A, tl)
    print(f"d={d} D={D} A={A}: |U/+-|={n} relations={len(ker)} max|r|_1={mx}: {name}  [{time.time()-t0:.1f}s]", flush=True)
