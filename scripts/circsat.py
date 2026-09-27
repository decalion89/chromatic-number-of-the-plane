"""Rational circular colourings by CP-SAT: characters phi = y/Q, y in Z^r, with every unit in the window.

Pure integer model:  Q*n_u + Q/k <= C_u . y <= Q*n_u + (k-1)Q/k   for every direction u,
y_i in [0, Q), n_u integer.  Finds every circular k-colouring of the module whose character has
denominator dividing Q (k must divide Q).  INFEASIBLE rules those out; it says nothing about other
denominators (use circgate.py's MILP for the continuous question).

usage: python3 circsat.py <module.json> <Q> [k] [seconds] [workers]"""
import sys, os, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path = [p for p in sys.path if os.path.abspath(p or ".") != HERE]   # no script here may shadow an installed package
from ortools.sat.python import cp_model
t0 = time.time()
src = open(os.path.join(HERE, "circgate.py")).read()
pre = src[:src.index("# optional: also require the multiples")]
ARGS = list(sys.argv)
sys.argv = ["circgate.py", ARGS[1], "1"]
exec(pre)                                   # (defines its own KC, TL: read ours afterwards)
path = ARGS[1]; Q = int(ARGS[2]); KC = int(ARGS[3]) if len(ARGS) > 3 else 5
TL = float(ARGS[4]) if len(ARGS) > 4 else 600; W = int(ARGS[5]) if len(ARGS) > 5 else 2
assert Q % KC == 0
reps, seenr = [], set()
for i in range(len(CU)):
    k_ = tuple(CU[i]); nk_ = tuple(-CU[i])
    if nk_ in seenr or k_ in seenr: continue
    seenr.add(k_); reps.append(i)
C = [[int(x) for x in CU[i]] for i in reps]
m = len(C)
mdl = cp_model.CpModel()
y = [mdl.NewIntVar(0, Q - 1, f"y{i}") for i in range(r)]
for j, cu in enumerate(C):
    lo = sum(min(0, c) for c in cu) * (Q - 1); hi = sum(max(0, c) for c in cu) * (Q - 1)
    n = mdl.NewIntVar(lo // Q - 1, hi // Q + 1, f"n{j}")
    expr = sum(c * yi for c, yi in zip(cu, y) if c)
    mdl.Add(expr - Q * n >= Q // KC)
    mdl.Add(expr - Q * n <= (KC - 1) * Q // KC)
sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = TL; sv.parameters.num_workers = W
st = sv.Solve(mdl)
print(f"{os.path.basename(path)}: {m} directions, rank {r}, Q = {Q}, k = {KC}: {sv.StatusName(st)}   [{time.time()-t0:.0f}s]", flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    yy = [sv.Value(v) for v in y]
    from fractions import Fraction as Fr
    print("  phi = " + " ".join(str(Fr(v, Q)) for v in yy))
