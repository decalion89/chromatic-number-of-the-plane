"""gen_kappa.py V.json N [A] [tlimit]: is there a character of ZU, U = G_N V (N = 5^k), with xi(U) in [A, 1-A]?
Units are written over the Q-basis of F (from the saved field) as integer 2n-vectors / D; relations by PARI matkerint;
MIP (SCIP via OR-tools): f_u in [A, 1-A], <r, f> integer for every relation r of a basis.  Float solver: a guide."""
import sys, json, math, subprocess, tempfile, os
from fractions import Fraction as Fr
from math import lcm
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from nf import Field
from ortools.linear_solver import pywraplp

def rotations(N):
    k = 0
    while 5 ** k < N: k += 1
    assert 5 ** k == N
    rho = (Fr(3, 5), Fr(4, 5)); out = []
    for j in range(-k, k + 1):
        g = (Fr(1), Fr(0))
        r = rho if j >= 0 else (rho[0], -rho[1])
        for _ in range(abs(j)):
            g = (g[0] * r[0] - g[1] * r[1], g[0] * r[1] + g[1] * r[0])
        out.append(g); out.append((-g[1], g[0]))          # g and i*g (one per +- pair)
    return out

def build(Fd, V, N):
    U = []
    for (X, Y) in V:
        xv, yv = Fd.vec(X), Fd.vec(Y)
        for (gr, gi) in rotations(N):
            U.append([gr * a - gi * b for a, b in zip(xv, yv)] + [gr * b + gi * a for a, b in zip(xv, yv)])
    D = 1
    for u in U:
        for x in u: D = lcm(D, x.denominator)
    return [[int(x * D) for x in u] for u in U], D

def relations(Ui):
    mat = "[" + ";".join(",".join(str(x) for x in u) for u in Ui) + "]"     # rows = units
    script = "A = " + mat + ";\nK = matkerint(A~);\nprint(#K);\nfor(j=1,#K, print(Vec(K[,j])));\nquit;\n"
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as fh:
        fh.write(script); fn = fh.name
    out = subprocess.run(["gp", "-q", "-s", "400000000", fn], capture_output=True, text=True, check=True).stdout.split("\n")
    os.unlink(fn)
    r = int(out[0])
    R = [[int(x) for x in out[1 + j].strip()[1:-1].split(",")] for j in range(r)]
    for rr in R:
        assert all(sum(rr[i] * Ui[i][c] for i in range(len(Ui))) == 0 for c in range(len(Ui[0])))
    return R

def solve(n, R, A, tlimit):
    s = pywraplp.Solver.CreateSolver("SCIP")
    f = [s.NumVar(float(A), float(1 - A), f"f{i}") for i in range(n)]
    for j, r in enumerate(R):
        lo = sum(c * (float(A) if c > 0 else float(1 - A)) for c in r)
        hi = sum(c * (float(1 - A) if c > 0 else float(A)) for c in r)
        z = s.IntVar(math.floor(lo) - 1, math.ceil(hi) + 1, f"z{j}")
        s.Add(sum(c * f[i] for i, c in enumerate(r) if c) == z)
    s.SetTimeLimit(int(tlimit * 1000))
    st = s.Solve()
    return {pywraplp.Solver.OPTIMAL: "FEASIBLE", pywraplp.Solver.FEASIBLE: "FEASIBLE",
            pywraplp.Solver.INFEASIBLE: "INFEASIBLE"}.get(st, f"UNKNOWN({st})")

if __name__ == "__main__":
    d = json.load(open(sys.argv[1])); N = int(sys.argv[2])
    A = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1, 3); tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
    Fd = Field([[Fr(c) for c in p] for p in d["polys"]])
    dec = lambda a: {tuple(int(x) for x in k.split(",")): Fr(v) for k, v in a.items()}
    V = [(dec(X), dec(Y)) for X, Y in d["V"]]
    Ui, D = build(Fd, V, N)
    R = relations(Ui)
    print(f"N = {N}: {len(Ui)} units (one per +- pair), D = {D}, {len(R)} relations, max |r| = {max(sum(abs(c) for c in r) for r in R)}: {solve(len(Ui), R, A, tl)}", flush=True)
