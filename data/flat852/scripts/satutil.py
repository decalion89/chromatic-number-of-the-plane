"""satutil.py -- pysat helpers (CaDiCaL 1.9.5 via pysat) with wall-clock limits via conflict-budget chunks,
and a DIMACS writer."""
import time
from pysat.solvers import Solver


def solve_limited(s, assumptions=(), tlimit=None, chunk=20000):
    """returns True/False/None (None = time out). CaDiCaL in pysat has no interrupt, so we run
    solve_limited with a conflict budget repeatedly (the solver keeps its state between calls)."""
    if tlimit is None:
        return s.solve(assumptions=list(assumptions))
    t0 = time.time()
    while True:
        s.conf_budget(chunk)
        r = s.solve_limited(assumptions=list(assumptions))
        if r is not None:
            return r
        if time.time() - t0 > tlimit:
            return None


def colouring_from_model(model, n, k=4):
    pos = set(l for l in model if l > 0)
    col = []
    for v in range(n):
        c = next((i for i in range(k) if (1 + k * v + i) in pos), None)
        col.append(c)
    return col


def check_colouring(col, E):
    return all(col[a] is not None and col[b] is not None and col[a] != col[b] for a, b in E)


def write_dimacs(path, nvars, clauses, comments=()):
    with open(path, "w") as f:
        for c in comments:
            f.write(f"c {c}\n")
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cl in clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")
