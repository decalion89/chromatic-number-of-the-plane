"""The circular gate in relation space.

Unknowns t_u = frac(phi(u)) in [1/k, 1 - 1/k], one per direction u (t_{-u} = 1 - t_u is implied).
Every integer relation sum_u a_u u = 0 among the unit vectors forces sum_u a_u t_u to be an integer.
Short relations (found by search: a_i u_i + a_j u_j + a_l u_l = 0 with small coefficients) give a MILP
with tiny integer ranges.  Using any set of relations is a RELAXATION: if it is infeasible, no
circular k-colouring exists along these units (F is empty).  If the relations have full rank n - r,
a solution t lifts to a character up to the finite index of the relation sublattice; the lift is
checked explicitly.

usage: python3 circrel.py <module.json> [k] [time limit] [max coefficient]"""
import sys, os, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
src = open(os.path.join(HERE, "circgate.py")).read()
pre = src[:src.index("# optional: also require the multiples")]
ARGS = list(sys.argv)
sys.argv = ["circgate.py", ARGS[1], "1"]
exec(pre)
KC = int(ARGS[2]) if len(ARGS) > 2 else 5; TLIM = float(ARGS[3]) if len(ARGS) > 3 else 600
AMAX = int(ARGS[4]) if len(ARGS) > 4 else 7
reps, seenr = [], set()
for i in range(len(CU)):
    k_ = tuple(CU[i]); nk_ = tuple(-CU[i])
    if nk_ in seenr or k_ in seenr: continue
    seenr.add(k_); reps.append(i)
D = CU[reps].astype(np.int64); n = len(D)
idx = {}
for j in range(n):
    idx[tuple(D[j])] = (j, 1); idx[tuple(-D[j])] = (j, -1)
FULL = os.environ.get("FULL", "1") == "1"
if FULL:
    # the whole relation lattice R = {a : a D = 0}, LLL-reduced: LLL on [I_n | W D]
    import flint
    Wt = 10 ** 12
    Mx = flint.fmpz_mat([[int(i == j) for j in range(n)] + [Wt * int(x) for x in D[i]] for i in range(n)])
    Lr = Mx.lll()
    rows = [[int(Lr[i, j]) for j in range(n + r)] for i in range(n)]
    kern = [row[:n] for row in rows if not any(row[n:])]
    rels = set()
    for row in kern:
        key = tuple((k2, v2) for k2, v2 in enumerate(row) if v2)
        neg = tuple((k2, -v2) for k2, v2 in key)
        rels.add(min(key, neg))
    print(f"  LLL kernel basis: {len(kern)} relations, max |coef| {max(max(abs(v) for v in row) for row in kern) if kern else 0}, max length {max(sum(1 for v in row if v) for row in kern) if kern else 0}   [{time.time()-t0:.0f}s]", flush=True)
    AMAX = 0
# search 3-term relations a*D_i + b*D_j = c*(sign) D_l
coefs = [c for c in range(-AMAX, AMAX + 1) if c]
rels = rels if FULL else set()
for i in range(n):
    for a in range(1, AMAX + 1):
        V = a * D[i] + np.outer(coefs, np.ones(r, dtype=np.int64))[:, None, :] * D[None, :, :]   # (coef, j, r)
        for ci, b in enumerate(coefs):
            for j in range(n):
                if j == i: continue
                v = V[ci, j]
                for c in range(1, AMAX + 1):
                    if np.any(v % c): continue
                    hit = idx.get(tuple(v // c))
                    if hit is None: continue
                    l, sg = hit
                    if l in (i, j): continue
                    # a D_i + b D_j - c*sg D_l = 0
                    rel = {i: a, j: b}; rel[l] = rel.get(l, 0) - c * sg
                    g = 0
                    for x in rel.values(): g = np.gcd(g, abs(x))
                    key = tuple(sorted((k2, v2 // g) for k2, v2 in rel.items() if v2))
                    neg = tuple(sorted((k2, -v2) for k2, v2 in key))
                    rels.add(min(key, neg))
# parallelogram relations: D_i + s_j D_j = D_k + s_l D_l (pairs of directions with the same sum)
from collections import defaultdict
sums = defaultdict(list) if not FULL else {}
for i in (range(n) if not FULL else []):
    for j in range(i + 1, n):
        for sj in (1, -1):
            sums[tuple(D[i] + sj * D[j])].append((i, j, sj))
for key_, lst in sums.items():
    if len(lst) < 2: continue
    i0, j0, s0 = lst[0]
    for (i1, j1, s1) in lst[1:]:
        rel = defaultdict(int)
        rel[i0] += 1; rel[j0] += s0; rel[i1] -= 1; rel[j1] -= s1
        key2 = tuple(sorted((k2, v2) for k2, v2 in rel.items() if v2))
        if not key2: continue
        neg = tuple(sorted((k2, -v2) for k2, v2 in key2))
        rels.add(min(key2, neg))
rels = sorted(rels)
Rm = np.zeros((len(rels), n), dtype=np.int64)
for q, rel in enumerate(rels):
    for k2, v2 in rel: Rm[q, k2] = v2
assert not (Rm @ D).any(), "a found relation is not a relation"
rk_rel = np.linalg.matrix_rank(Rm.astype(float)) if len(rels) else 0
print(f"{os.path.basename(ARGS[1])}: {n} directions, rank {r}; {len(rels)} short relations (|coef| <= {AMAX}), rank {rk_rel} of {n - r} needed   [{time.time()-t0:.0f}s]", flush=True)
# MILP: t in [1/k + s, 1 - 1/k - s]; for each relation: sum a_u t_u - z = 0, z integer in its range
lo, hi = 1.0 / KC, 1.0 - 1.0 / KC
nrel = len(rels); nv = n + nrel + 1
A = np.zeros((nrel, nv)); A[:, :n] = Rm; A[np.arange(nrel), n + np.arange(nrel)] = -1.0
zlo = np.floor(np.where(Rm > 0, Rm * lo, Rm * hi).sum(axis=1)) ; zhi = np.ceil(np.where(Rm > 0, Rm * hi, Rm * lo).sum(axis=1))
Bt1 = np.zeros((n, nv)); Bt1[np.arange(n), np.arange(n)] = 1.0; Bt1[:, -1] = -1.0      # t - s >= lo
Bt2 = np.zeros((n, nv)); Bt2[np.arange(n), np.arange(n)] = 1.0; Bt2[:, -1] = 1.0       # t + s <= hi
cons = [LinearConstraint(A, lb=0, ub=0), LinearConstraint(Bt1, lb=lo, ub=np.inf), LinearConstraint(Bt2, lb=-np.inf, ub=hi)]
lb = np.concatenate([np.full(n, lo), zlo, [0.0]]); ub = np.concatenate([np.full(n, hi), zhi, [0.5]])
integ = np.concatenate([np.zeros(n), np.ones(nrel), [0]])
cobj = np.zeros(nv); cobj[-1] = -1.0
SOLVER = os.environ.get("SOLVER", "HIGHS")
if SOLVER == "HIGHS":
    res = milp(cobj, constraints=cons, integrality=integ, bounds=Bounds(lb, ub), options={"time_limit": TLIM})
else:
    # an independent branch and bound (SCIP or CBC through OR-Tools) on the same model, as a cross-check
    from types import SimpleNamespace
    from ortools.linear_solver import pywraplp
    so = pywraplp.Solver.CreateSolver(SOLVER)
    so.SetTimeLimit(int(TLIM * 1000))
    tv = [so.NumVar(lo, hi, f"t{j}") for j in range(n)]
    zv = [so.IntVar(float(zlo[q]), float(zhi[q]), f"z{q}") for q in range(nrel)]
    sv = so.NumVar(0.0, 0.5, "s")
    for q, rel in enumerate(rels):
        so.Add(sum(int(v2) * tv[k2] for k2, v2 in rel) == zv[q])
    for j in range(n):
        so.Add(tv[j] - sv >= lo); so.Add(tv[j] + sv <= hi)
    so.Maximize(sv)
    st = so.Solve()
    code = {pywraplp.Solver.OPTIMAL: 0, pywraplp.Solver.FEASIBLE: 1, pywraplp.Solver.INFEASIBLE: 2}.get(st, 4)
    xs = np.array([v.solution_value() for v in tv] + [0.0] * nrel + [sv.solution_value()]) if code in (0, 1) else None
    res = SimpleNamespace(status=code, message=f"{SOLVER} status {st}", x=xs)
print(f"  relation-space MILP: status {res.status}: {res.message}   [{time.time()-t0:.0f}s]", flush=True)
if res.status == 2:
    print(f"  ==> INFEASIBLE: no circular {KC}-colouring along these units (the relations used are genuine, so this is a relaxation)")
elif res.x is not None:
    t = res.x[:n]
    print(f"  solution with slack {res.x[-1]:.6f}")
    viol = np.abs(((Rm @ t) + 0.5) % 1.0 - 0.5).max() if nrel else 0.0
    full = FULL and rk_rel == n - r
    print(f"  relations satisfied up to {viol:.2e}; relation lattice {'COMPLETE: t is a genuine character' if full else 'incomplete: t may not lift'}")
