"""per_ilp.py REL.json [REL2.json ...]: is there an integer homomorphism p on the relation lattice (p_j = p(r_j)) with
p(rho) strictly inside (P/4 - 3N/4, 3P/4 - N/4) for every listed relation rho (P, N: sums of the positive and negative
coefficients)?  This is the base-point period system of Lemma F17 when the closed walks of the relations and enough
squares lie in H.  Infeasible -> prints a small infeasible core (CP-SAT assumptions)."""
import json, sys
from fractions import Fraction as Fr
from math import floor, ceil
from ortools.sat.python import cp_model
rels = []
for f in sys.argv[1:]:
    rels += json.load(open(f))
seen = set(); R = []
for r in rels:
    k = tuple(r["rho"]); km = tuple(-x for x in k)
    if k in seen or km in seen:
        continue
    seen.add(k); R.append(r)
n = len(R[0]["a"])
md = cp_model.CpModel()
p = [md.NewIntVar(-1000, 1000, f"p{j}") for j in range(n)]
lits = []
for i, r in enumerate(R):
    P = sum(x for x in r["rho"] if x > 0); N = -sum(x for x in r["rho"] if x < 0)
    A = Fr(P - 3 * N, 4); B = Fr(3 * P - N, 4)
    lo, hi = floor(A) + 1, ceil(B) - 1
    b = md.NewBoolVar(f"c{i}"); lits.append(b)
    e = sum(a * p[j] for j, a in enumerate(r["a"]) if a)
    md.Add(e >= lo).OnlyEnforceIf(b); md.Add(e <= hi).OnlyEnforceIf(b)
md.AddAssumptions(lits)
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 600; s.parameters.num_workers = 1
st = s.Solve(md)
print(len(R), "relations:", s.StatusName(st))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    print("p =", [s.Value(x) for x in p])
elif st == cp_model.INFEASIBLE:
    core = s.SufficientAssumptionsForInfeasibility()
    idx = sorted(int(md.Proto().variables[c].name[1:]) for c in core)
    print("core of", len(idx), "relations; l1 norms:", sorted(R[i]["l1"] for i in idx))
    json.dump([R[i] for i in idx], open("core.json", "w"))
