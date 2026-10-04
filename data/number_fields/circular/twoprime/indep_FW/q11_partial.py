"""Partial lower-bound checks for the 76-vertex graph: UNSAT at fractions just below 16/5 (own encoding, pysat)."""
import os
import json, time, sys
from pysat.solvers import Solver
g = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "quadratic_planes", "q11.json")))
E = [tuple(e) for e in g["edges"]]; n = len(g["points"])
def hom(p, q):
    var = lambda v, a: v * p + a + 1
    s = Solver(name="cadical153")
    for v in range(n):
        s.add_clause([var(v, a) for a in range(p)])
    bad = [d for d in range(p) if not (q <= d <= p - q)]
    for (u, v) in E:
        for a in range(p):
            for d in bad:
                s.add_clause([-var(u, a), -var(v, (a + d) % p)])
    s.add_clause([var(E[0][0], 0)])
    t = time.time(); ok = s.solve(); dt = time.time() - t; s.delete(); return ok, dt
for p, q in [(int(a), int(b)) for a, b in (x.split("/") for x in sys.argv[1:])]:
    ok, dt = hom(p, q)
    print(f"K_{p}/{q} ({p/q:.5f}): {'SAT' if ok else 'UNSAT'} [{dt:.1f}s]", flush=True)
