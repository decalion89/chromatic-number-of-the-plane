"""kappaD.py d D [A] [tlimit]: Lemma W test on U_D = unit vectors ((a + b sqrt d)/D, (c + e sqrt d)/D) of Q(sqrt d)^2.
Characters of Z U_D are restrictions of theta in T^4 acting on (1/D)Z^4: xi(u) = <theta, (a, b, c, e)>.
Question: is there theta with <theta, v_u> in [A, 1-A] (mod 1) for every u?  For A = 1/3 this holds iff
Cay(Z U_D, U_D) is 3-colourable (Lemma W).  Numerical MILP (SCIP via OR-tools); certify separately."""
import sys, time
from math import isqrt
from fractions import Fraction as Fr
from ortools.linear_solver import pywraplp

def units(d, D):
    out = set()
    D2 = D * D
    bmax = isqrt(D2 // d)
    for b in range(-bmax, bmax + 1):
        for e in range(-bmax, bmax + 1):
            rest = D2 - d * (b * b + e * e)
            if rest < 0: continue
            if b == 0 and e == 0:
                for a in range(-D, D + 1):
                    c2 = rest - a * a
                    if c2 < 0: continue
                    c = isqrt(c2)
                    if c * c == c2:
                        out.add((a, 0, c, 0)); out.add((a, 0, -c, 0))
                continue
            # (a, c) = lam * (-e, b), lam^2 (b^2 + e^2) = rest
            n2 = b * b + e * e
            if rest % 1: continue
            q = Fr(rest, n2)
            num, den = q.numerator, q.denominator
            rn, rd = isqrt(num), isqrt(den)
            if rn * rn != num or rd * rd != den: continue
            lam = Fr(rn, rd)
            for s in (1, -1):
                a, c = s * lam * (-e), s * lam * b
                if a.denominator == 1 and c.denominator == 1:
                    out.add((int(a), b, int(c), e))
    for (a, b, c, e) in out:
        assert a * a + d * b * b + c * c + d * e * e == D2 and a * b + c * e == 0
    return sorted(out)

def halve(U):
    seen = set(); H = []
    for u in U:
        if tuple(-x for x in u) in seen: continue
        seen.add(u); H.append(u)
    return H

def solve(V, A, tlimit=600, extra=None):
    s = pywraplp.Solver.CreateSolver("SCIP")
    th = [s.NumVar(0, 1, f"t{j}") for j in range(4)]
    for i, v in enumerate(V):
        lo = sum(min(0, c) for c in v) - 1; hi = sum(max(0, c) for c in v) + 1
        k = s.IntVar(lo, hi, f"k{i}")
        expr = sum(c * th[j] for j, c in enumerate(v) if c) - k
        s.Add(expr >= A); s.Add(expr <= 1 - A)
    s.SetTimeLimit(int(tlimit * 1000))
    st = s.Solve()
    name = {pywraplp.Solver.OPTIMAL: "FEASIBLE", pywraplp.Solver.FEASIBLE: "FEASIBLE",
            pywraplp.Solver.INFEASIBLE: "INFEASIBLE"}.get(st, f"UNKNOWN(status {st})")
    return name, ([t.solution_value() for t in th] if name == "FEASIBLE" else None)

if __name__ == "__main__":
    d, D = int(sys.argv[1]), int(sys.argv[2])
    A = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1, 3)
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    t0 = time.time()
    U = halve(units(d, D))
    name, th = solve(U, float(A), tl)
    print(f"d={d} D={D} A={A}: |U_D/+-|={len(U)}: {name}  [{time.time()-t0:.1f}s]", flush=True)
    if th:
        vals = [sum(c * t for c, t in zip(v, th)) % 1 for v in U]
        print("   theta =", [round(x, 6) for x in th], " min dist to 0:", round(min(min(x, 1 - x) for x in vals), 6))
