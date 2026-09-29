"""phase_solve.py PTS.npy STATE.json TLIMIT -- CaDiCaL (pysat) 4-colouring with phases from a partial colouring."""
import sys, json, time
import numpy as np
from flat import *
from satutil import solve_limited
from pysat.solvers import Solver

D, CD, U = directions()
P = np.load(sys.argv[1]); st = json.load(open(sys.argv[2])); TL = float(sys.argv[3])
tri = st["triangle"]; col = list(st.get("col_partial") or st.get("col"))
n = len(P)
col += [-1] * (n - len(col))
E, J = build_edges(P, U)
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(int(b)); adj[b].add(int(a))
for v in range(n):            # greedy phase for uncoloured vertices
    if col[v] < 0:
        used = [sum(1 for w in adj[v] if col[w] == c) for c in range(4)]
        col[v] = int(np.argmin(used))
for i, v in enumerate(tri):
    col[v] = i
conf = sum(1 for a, b in E if col[a] == col[b])
s = Solver(name="cadical195", bootstrap_with=cnf_clauses(n, E.tolist(), tri))
s.set_phases([1 + 4 * v + col[v] for v in range(n)] + [-(1 + 4 * v + c) for v in range(n) for c in range(4) if c != col[v]])
print(f"n={n} m={len(E)}; initial phase colouring has {conf} conflicts", flush=True)
t = time.time()
r = solve_limited(s, tlimit=TL, chunk=50000)
dt = time.time() - t
if r:
    m = set(l for l in s.get_model() if l > 0)
    c = [next(k for k in range(4) if 1 + 4 * v + k in m) for v in range(n)]
    ok = all(c[a] != c[b] for a, b in E)
    json.dump({"triangle": tri, "col": c}, open(sys.argv[1].replace("_pts.npy", "_colour.json"), "w"))
    print(f"SAT in {dt:.0f}s; colouring verified: {ok}", flush=True)
else:
    print(f"{'UNSAT' if r is False else 'UNKNOWN'} after {dt:.0f}s", flush=True)
