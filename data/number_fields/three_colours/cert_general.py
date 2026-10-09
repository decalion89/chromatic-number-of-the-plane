"""cert_general.py V.json N OUT.json: exact certificate that no character of ZU, U = G_N V (one unit per +- pair),
maps U into [1/3, 2/3], for a multiquadratic field (the branch and bound of certify_w2.py, any dimension).
Relations: PARI matkerint (checked exactly).  Floats (GLOP) only guide the search; every leaf is checked exactly."""
import sys, json, time, math
from fractions import Fraction as Fr
from ortools.linear_solver import pywraplp
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from nf import Field
from gen_kappa import build, relations

LO, HI = Fr(1, 3), Fr(2, 3)

def box_range(r):
    lo = sum(c * (LO if c > 0 else HI) for c in r); hi = sum(c * (HI if c > 0 else LO) for c in r)
    return math.ceil(lo), math.floor(hi)

def exact_farkas_ok(rows, rhs, y):
    n = len(rows[0]); c = [sum(y[k] * rows[k][i] for k in range(len(rows))) for i in range(n)]
    mx = sum(ci * (HI if ci > 0 else LO) for ci in c); mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
    t = sum(y[k] * rhs[k] for k in range(len(rows)))
    return t > mx or t < mn

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
        if s.Solve() != pywraplp.Solver.OPTIMAL or s.Objective().Value() <= 1e-9:
            continue
        yv = [v.solution_value() for v in y]
        for lim in (10 ** 3, 10 ** 5, 10 ** 7, 10 ** 9):
            yr = [Fr(v).limit_denominator(lim) for v in yv]
            if exact_farkas_ok(rows, rhs, yr):
                return yr
    return None

def lp_feasible(rows, rhs, n):
    s = pywraplp.Solver.CreateSolver("GLOP")
    f = [s.NumVar(float(LO), float(HI), f"f{i}") for i in range(n)]
    for r, z in zip(rows, rhs):
        s.Add(sum(c * f[i] for i, c in enumerate(r) if c) == z)
    return s.Solve() in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE)

def lp_range(rows, rhs, n, r):
    out = []
    for sense in (1, -1):
        s = pywraplp.Solver.CreateSolver("GLOP")
        f = [s.NumVar(float(LO), float(HI), f"f{i}") for i in range(n)]
        for rr, z in zip(rows, rhs):
            s.Add(sum(c * f[i] for i, c in enumerate(rr) if c) == z)
        s.Minimize(sense * sum(c * f[i] for i, c in enumerate(r) if c))
        if s.Solve() != pywraplp.Solver.OPTIMAL:
            return None
        out.append(sense * s.Objective().Value())
    return out[0], out[1]

def main():
    d = json.load(open(sys.argv[1])); N = int(sys.argv[2])
    Fd = Field([[Fr(c) for c in p] for p in d["polys"]])
    dec = lambda a: {tuple(int(x) for x in k.split(",")): Fr(v) for k, v in a.items()}
    V = [(dec(X), dec(Y)) for X, Y in d["V"]]
    U, D = build(Fd, V, N); n = len(U)
    R = relations(U)
    R.sort(key=lambda r: sum(abs(c) for c in r))
    stats = {"nodes": 0, "leaves": 0, "last": time.time()}; t0 = time.time()

    def node(fixed):
        stats["nodes"] += 1
        if time.time() - stats["last"] > 30:
            stats["last"] = time.time()
            print(f"  nodes {stats['nodes']} leaves {stats['leaves']} depth {len(fixed)} [{time.time()-t0:.0f}s]", flush=True)
        rows = [R[j] for j, _ in fixed]; rhs = [z for _, z in fixed]
        if rows and not lp_feasible(rows, rhs, n):
            y = farkas(rows, rhs)
            if y is None:
                raise RuntimeError(f"no exact Farkas vector at depth {len(fixed)}")
            stats["leaves"] += 1
            return {"leaf": [str(v) for v in y]}
        done = {j for j, _ in fixed}
        free = [j for j in range(len(R)) if j not in done]
        if not free:
            raise RuntimeError("feasible leaf: a candidate character exists")
        best = None
        for j in free:
            rg = lp_range(rows, rhs, n, R[j])
            cnt = 0 if rg is None else max(0, math.floor(rg[1] + 1e-9) - math.ceil(rg[0] - 1e-9) + 1)
            if best is None or cnt < best[0]:
                best = (cnt, j)
            if cnt == 0:
                break
        j = best[1]
        lo, hi = box_range(R[j])
        kids = {}
        for z in range(lo, hi + 1):
            kids[str(z)] = node(fixed + [(j, z)])
        return {"branch": j, "lo": lo, "hi": hi, "kids": kids}

    tree = node([])
    gens = [int(-Fr(p[0])) for p in d["polys"]]
    json.dump({"gens": gens, "basis": [list(b) for b in Fd.basis], "N": N, "D": D, "units": U, "relations": R,
               "tree": tree}, open(sys.argv[3], "w"))
    print(f"certificate written: {n} units, {len(R)} relations, {stats['nodes']} nodes, {stats['leaves']} leaves  [{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
