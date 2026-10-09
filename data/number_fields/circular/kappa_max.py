"""kappa(U) = max over characters xi of ZU of min_u ||xi(u)||, U = G_N V, for Q(sqrt d) with V = {1, u_n, conj u_n
(n in NS)}: MIP (SCIP via OR-tools): maximise t with f_u in [t, 1 - t] (one f per +- pair) and <r, f> integer for a
basis of the integer relations.  chi_c(Cay(ZU, U)) = 1/kappa when kappa > 1/4 (Theorem W+)."""
import sys, math
from fractions import Fraction as Fr
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "three_colours"))
from nf import Field, unit_from, lconj
from gen_kappa import build, relations
from ortools.linear_solver import pywraplp
d = int(sys.argv[1]); N = int(sys.argv[2]); ns = [int(x) for x in sys.argv[3].split(",")]; tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600
Fd = Field([[-d, 0, 1]]); sq = Fd.gen(0)
V = [(Fd.one(), {})]
for n in ns:
    u = unit_from(Fd, Fd.scal(Fr(n, d), sq)); V += [u, lconj(Fd, u)]
U, D = build(Fd, V, N); R = relations(U); n = len(U)
s = pywraplp.Solver.CreateSolver("SCIP")
t = s.NumVar(0.25, 0.5, "t")
f = [s.NumVar(0, 1, f"f{i}") for i in range(n)]
for i in range(n):
    s.Add(f[i] >= t); s.Add(f[i] <= 1 - t)
for j, r in enumerate(R):
    lo = sum(c * (0 if c > 0 else 1) for c in r); hi = sum(c * (1 if c > 0 else 0) for c in r)
    z = s.IntVar(math.floor(lo) - 1, math.ceil(hi) + 1, f"z{j}")
    s.Add(sum(c * f[i] for i, c in enumerate(r) if c) == z)
s.Maximize(t); s.SetTimeLimit(int(tl * 1000))
st = s.Solve()
print(f"d = {d}, N = {N}, n in {ns}: |U| = {n}: status {st}, kappa ~ {t.solution_value():.6f} (2/7 = {2/7:.6f}), bound {s.Objective().BestBound():.6f}", flush=True)
