"""cert_open_units.py IN.json OUT.json P Q [maxnodes] [eps] [ncand] [prefix]: exact certificate that no character of
Z U maps every vector of a given finite set U of unit vectors into the OPEN interval (Q/P, 1 - Q/P), 2Q < P <= 4Q.
IN.json: {"d": d, "D": D, "units": [[...], ...]} (one vector per +- pair, integer coordinates over D: 4 of them,
(x0 + x1 sqrt d)/D and (y0 + y1 sqrt d)/D, for Q(sqrt d); 8, over (1, sqrt3, sqrt11, sqrt33), for d = "3,11").
Relations: PARI matkerint (an LLL-reduced basis), sorted by l1 norm.  Branch on the integers z_j = <r_j, f> strictly
inside the open range; at a leaf an exact Farkas vector y with <y, z> >= max or <= min of <yR, f> over the CLOSED box
proves that no f in the open box solves the fixed equations (if yR = 0, <y, z> != 0).  Same method and output format
as ../cert_open.py ("N": 0); the units are checked by the checkers (../check_open.py, ../check_open_indep.py), not
here.  With a prefix (a JSON list of [relation, value] pairs) it writes the certificate of that subtree only."""
import sys, json, time, math
from fractions import Fraction as Fr
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..", "three_colours"))
from gen_kappa import relations
from ortools.linear_solver import pywraplp
sys.setrecursionlimit(100000)

inp, out, P, Q = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
maxnodes = int(sys.argv[5]) if len(sys.argv) > 5 else 10 ** 7
EPSARG = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-7
NCAND = int(sys.argv[7]) if len(sys.argv) > 7 else 10 ** 9
PREFIX = [tuple(x) for x in json.loads(sys.argv[8])] if len(sys.argv) > 8 else []
C = json.load(open(inp)); d, D, U = C["d"], C["D"], C["units"]
assert 2 * Q < P <= 4 * Q
LO, HI = Fr(Q, P), 1 - Fr(Q, P)
EPS = EPSARG        # float guide: closed box [LO + EPS, HI - EPS]
n = len(U)

def open_range(r):
    mn = sum(c * (LO if c > 0 else HI) for c in r); mx = sum(c * (HI if c > 0 else LO) for c in r)
    return math.floor(mn) + 1, math.ceil(mx) - 1

def leaf_ok(rows, rhs, y):
    c = [sum(y[k] * rows[k][i] for k in range(len(rows))) for i in range(n)]
    t = sum(y[k] * rhs[k] for k in range(len(rows)))
    if all(ci == 0 for ci in c):
        return t != 0
    mx = sum(ci * (HI if ci > 0 else LO) for ci in c); mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
    return t >= mx or t <= mn

def farkas(rows, rhs):
    m = len(rows)
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

def mk_lp(rows, rhs):
    s = pywraplp.Solver.CreateSolver("GLOP")
    f = [s.NumVar(float(LO) + EPS, float(HI) - EPS, f"f{i}") for i in range(n)]
    for r, z in zip(rows, rhs):
        s.Add(sum(c * f[i] for i, c in enumerate(r) if c) == z)
    return s, f

def lp_feasible(rows, rhs):
    s, f = mk_lp(rows, rhs)
    return s.Solve() in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE)

def lp_range(model, r):
    s, f = model
    res = []
    for sense in (1, -1):
        obj = s.Objective(); obj.Clear()
        for i, c in enumerate(r):
            if c: obj.SetCoefficient(f[i], sense * c)
        obj.SetMinimization()
        if s.Solve() != pywraplp.Solver.OPTIMAL:
            return None
        res.append(sense * s.Objective().Value())
    return res[0], res[1]

t0 = time.time()
R = relations(U); R.sort(key=lambda r: sum(abs(c) for c in r))
print(f"d={d} D={D}: {n} units, {len(R)} relations, open interval ({Q}/{P}, {P - Q}/{P}) [{time.time()-t0:.1f}s]", flush=True)
stats = {"nodes": 0, "leaves": 0, "last": time.time()}

def node(fixed):
    stats["nodes"] += 1
    if stats["nodes"] > maxnodes:
        raise RuntimeError("node limit")
    if time.time() - stats["last"] > 60:
        stats["last"] = time.time()
        print(f"  nodes {stats['nodes']} leaves {stats['leaves']} depth {len(fixed)} [{time.time()-t0:.0f}s]", flush=True)
    rows = [R[j] for j, _ in fixed]; rhs = [z for _, z in fixed]
    if rows and not lp_feasible(rows, rhs):
        y = farkas(rows, rhs)
        if y is not None:
            stats["leaves"] += 1
            return {"leaf": [str(v) for v in y]}
    done = {j for j, _ in fixed}
    free = [j for j in range(len(R)) if j not in done]
    if not free:
        raise RuntimeError(f"feasible leaf at {fixed}: a character into the open interval may exist")
    best = None
    model = mk_lp(rows, rhs)
    for j in free[:NCAND]:
        rg = lp_range(model, R[j])
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

tree = node(PREFIX)
if PREFIX:
    json.dump({"d": d, "D": D, "P": P, "Q": Q, "prefix": PREFIX, "subtree": tree, "nrel": len(R)}, open(out, "w"))
else:
    json.dump({"d": d, "N": 0, "D": D, "P": P, "Q": Q, "units": U, "relations": R, "tree": tree}, open(out, "w"))
print(f"certificate written: {n} units, {len(R)} relations, {stats['nodes']} nodes, {stats['leaves']} leaves [{time.time()-t0:.0f}s]", flush=True)
