"""per_cegar.py REL.json ITER OUT.json: counterexample-guided search for a finite set S of relations whose base-point
period system (Lemma F17) is infeasible.  S starts as the relations in REL.json; each round solves for an integer
homomorphism p (p_j = p(r_j)) inside every open range of S, then adds the shortest relation w (l1 norm, CP-SAT) with
p(w) outside its open range.  Prints the lengths of the added relations and the LP margin t* of each p."""
import json, sys, time, gzip, os
from fractions import Fraction as Fr
from math import floor, ceil
from ortools.sat.python import cp_model
import numpy as np
from scipy.optimize import linprog

C = json.load(gzip.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "at_four", "cert311_open_4.json.gz"), "rt")); R = C["relations"]; n, m = len(R), len(R[0])
Rm = np.array(R, dtype=float)
S = json.load(open(sys.argv[1])); ITER = int(sys.argv[2]); out = sys.argv[3]
t0 = time.time()


def rng(rho):
    P = sum(x for x in rho if x > 0); N = -sum(x for x in rho if x < 0)
    return floor(Fr(P - 3 * N, 4)) + 1, ceil(Fr(3 * P - N, 4)) - 1


def solve_p(S):
    md = cp_model.CpModel()
    p = [md.NewIntVar(-200, 200, f"p{j}") for j in range(n)]
    for r in S:
        lo, hi = rng(r["rho"])
        e = sum(a * p[j] for j, a in enumerate(r["a"]) if a)
        md.Add(e >= lo); md.Add(e <= hi)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 1200; s.parameters.num_workers = 1
    st = s.Solve(md)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return [s.Value(x) for x in p]
    return s.StatusName(st)


def margin(p):
    """max t with R f = p, 1/4 + t <= f_u <= 3/4 - t (LP, float)"""
    c = np.zeros(m + 1); c[-1] = -1
    Aeq = np.hstack([Rm, np.zeros((n, 1))]); beq = np.array(p, dtype=float)
    Aub = np.vstack([np.hstack([-np.eye(m), np.ones((m, 1))]), np.hstack([np.eye(m), np.ones((m, 1))])])
    bub = np.concatenate([-0.25 * np.ones(m), 0.75 * np.ones(m)])
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(None, None)] * (m + 1), method="highs")
    return res.x[-1] if res.status == 0 else None


def refute(p, K=40):
    md = cp_model.CpModel()
    y = [md.NewIntVar(-K, K, f"y{j}") for j in range(n)]
    wp = [md.NewIntVar(0, 10 * K, f"wp{u}") for u in range(m)]
    wn = [md.NewIntVar(0, 10 * K, f"wn{u}") for u in range(m)]
    for u in range(m):
        md.Add(sum(R[j][u] * y[j] for j in range(n) if R[j][u]) == wp[u] - wn[u])
    P = sum(wp); N = sum(wn)
    md.Add(P + N >= 1)
    md.Add(4 * sum(p[j] * y[j] for j in range(n) if p[j]) >= 3 * P - N)
    md.Minimize(P + N)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 600; s.parameters.num_workers = 1
    st = s.Solve(md)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, s.StatusName(st)
    a = [s.Value(v) for v in y]
    rho = [sum(a[j] * R[j][u] for j in range(n)) for u in range(m)]
    return {"a": a, "rho": rho, "l1": sum(abs(x) for x in rho)}, s.StatusName(st)


added = []
for it in range(ITER):
    p = solve_p(S)
    if isinstance(p, str):
        print(f"round {it}: {p} with {len(S)} relations [{time.time()-t0:.0f}s]", flush=True)
        break
    t = margin(p)
    w, st = refute(p)
    if w is None:
        print(f"round {it}: p feasible, no refutation found ({st}); t* = {t}", flush=True)
        break
    S.append(w); added.append(w)
    print(f"round {it}: t* = {t:.5f}, added relation of l1 {w['l1']} ({st}) [{time.time()-t0:.0f}s]", flush=True)
    json.dump(added, open(out, "w"))
