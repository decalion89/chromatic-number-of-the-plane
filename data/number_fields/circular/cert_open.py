"""cert_open.py OUT.json d N n1,n2,... P Q : exact certificate that no character of ZU (U = G_N V, V = {1, u_n, conj u_n}
over Q(sqrt d)) maps every u into the OPEN interval (Q/P, 1 - Q/P).  Then kappa(U) <= Q/P, and by Theorem W+ the
Cayley graph has no homomorphism to any K_{p/q} with p/q < P/Q, i.e. chi_c >= P/Q.
Branch on the integers z_j = <r_j, f> strictly inside the open range; at a leaf an exact Farkas vector y with
<y, z> >= max or <= min of <yR, f> over the CLOSED box (f in [Q/P, 1 - Q/P]) proves that no f in the open box solves
the fixed equations (if yR = 0 we need <y, z> != 0)."""
import sys, json, time, math
from fractions import Fraction as Fr
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "three_colours"))
from nf import Field, unit_from, lconj
from gen_kappa import build, relations
from ortools.linear_solver import pywraplp

out, d, N, ns, P, Q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), [int(x) for x in sys.argv[4].split(",")], int(sys.argv[5]), int(sys.argv[6])
LO, HI = Fr(Q, P), 1 - Fr(Q, P)
EPS = 1e-7          # float guide: closed box [LO + EPS, HI - EPS]

def open_range(r):
    """integers strictly between the min and max of <r, f> over the closed box"""
    mn = sum(c * (LO if c > 0 else HI) for c in r); mx = sum(c * (HI if c > 0 else LO) for c in r)
    lo = math.floor(mn) + 1; hi = math.ceil(mx) - 1
    return lo, hi

def leaf_ok(rows, rhs, y):
    n = len(rows[0]); c = [sum(y[k] * rows[k][i] for k in range(len(rows))) for i in range(n)]
    t = sum(y[k] * rhs[k] for k in range(len(rows)))
    if all(ci == 0 for ci in c):
        return t != 0
    mx = sum(ci * (HI if ci > 0 else LO) for ci in c); mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
    return t >= mx or t <= mn

def farkas(rows, rhs):
    n = len(rows[0]); m = len(rows)
    for sgn in (1, -1):
        s = pywraplp.Solver.CreateSolver("GLOP")
        y = [s.NumVar(-1, 1, f"y{k}") for k in range(m)]
        tt = [s.NumVar(-s.infinity(), s.infinity(), f"t{i}") for i in range(n)]
        for i in range(n):
            ci = sum(rows[k][i] * y[k] for k in range(m) if rows[k][i])
            s.Add(tt[i] >= sgn * ci * float(LO)); s.Add(tt[i] >= sgn * ci * float(HI))
        s.Maximize(sgn * sum(rhs[k] * y[k] for k in range(m)) - sum(tt))
        if s.Solve() != pywraplp.Solver.OPTIMAL or s.Objective().Value() < -1e-9:
            continue
        yv = [v.solution_value() for v in y]
        for lim in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 7, 10 ** 9):
            yr = [Fr(v).limit_denominator(lim) for v in yv]
            if leaf_ok(rows, rhs, yr):
                return yr
    return None

def lp_feasible(rows, rhs, n):
    s = pywraplp.Solver.CreateSolver("GLOP")
    f = [s.NumVar(float(LO) + EPS, float(HI) - EPS, f"f{i}") for i in range(n)]
    for r, z in zip(rows, rhs):
        s.Add(sum(c * f[i] for i, c in enumerate(r) if c) == z)
    return s.Solve() in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE)

def lp_range(rows, rhs, n, r):
    res = []
    for sense in (1, -1):
        s = pywraplp.Solver.CreateSolver("GLOP")
        f = [s.NumVar(float(LO) + EPS, float(HI) - EPS, f"f{i}") for i in range(n)]
        for rr, z in zip(rows, rhs):
            s.Add(sum(c * f[i] for i, c in enumerate(rr) if c) == z)
        s.Minimize(sense * sum(c * f[i] for i, c in enumerate(r) if c))
        if s.Solve() != pywraplp.Solver.OPTIMAL:
            return None
        res.append(sense * s.Objective().Value())
    return res[0], res[1]

Fd = Field([[-d, 0, 1]]); sq = Fd.gen(0)
V = [(Fd.one(), {})]
for nn in ns:
    u = unit_from(Fd, Fd.scal(Fr(nn, d), sq)); V += [u, lconj(Fd, u)]
U, D = build(Fd, V, N); n = len(U)
R = relations(U); R.sort(key=lambda r: sum(abs(c) for c in r))
stats = {"nodes": 0, "leaves": 0, "last": time.time()}; t0 = time.time()

def node(fixed):
    stats["nodes"] += 1
    if time.time() - stats["last"] > 30:
        stats["last"] = time.time()
        print(f"  nodes {stats['nodes']} leaves {stats['leaves']} depth {len(fixed)} [{time.time()-t0:.0f}s]", flush=True)
    rows = [R[j] for j, _ in fixed]; rhs = [z for _, z in fixed]
    if rows and not lp_feasible(rows, rhs, n):
        y = farkas(rows, rhs)
        if y is not None:
            stats["leaves"] += 1
            return {"leaf": [str(v) for v in y]}
    done = {j for j, _ in fixed}
    free = [j for j in range(len(R)) if j not in done]
    if not free:
        raise RuntimeError("feasible leaf: a character into the open interval may exist")
    best = None
    for j in free:
        rg = lp_range(rows, rhs, n, R[j])
        cnt = 0 if rg is None else max(0, math.floor(rg[1] + 1e-9) - math.ceil(rg[0] - 1e-9) + 1)
        if best is None or cnt < best[0]:
            best = (cnt, j)
        if cnt == 0: break
    j = best[1]
    lo, hi = open_range(R[j])
    kids = {}
    for z in range(lo, hi + 1):
        kids[str(z)] = node(fixed + [(j, z)])
    return {"branch": j, "lo": lo, "hi": hi, "kids": kids}

tree = node([])
json.dump({"d": d, "N": N, "D": D, "P": P, "Q": Q, "units": U, "relations": R, "tree": tree}, open(out, "w"))
print(f"certificate written: {n} units, {len(R)} relations, {stats['nodes']} nodes, {stats['leaves']} leaves [{time.time()-t0:.0f}s]")
