"""lns.py -- large-neighbourhood repair of a 4-colouring.
Given a graph, a colouring of the old vertices and uncoloured new vertices: free the vertices within BFS
distance r of the new ones, fix all others to their old colours, solve the local problem (CaDiCaL via pysat);
if UNSAT, increase r. Returns a full proper colouring or None."""
import time
import numpy as np
from pysat.solvers import Solver
from satutil import solve_limited


def lns_repair(n, adj, col, new, tri, rmax=6, tl=300, log=print):
    col = list(col) + [-1] * (n - len(col))
    for i, v in enumerate(tri):
        col[v] = i
    dist = {v: 0 for v in new}
    frontier = list(new)
    r = 0
    while True:
        # region: dist <= r
        R = [v for v, d in dist.items() if d <= r and v not in tri]
        Rs = set(R)
        idx = {v: i for i, v in enumerate(R)}
        var = lambda v, c: 1 + 4 * idx[v] + c
        s = Solver(name="cadical195")
        for v in R:
            s.add_clause([var(v, c) for c in range(4)])
            for w in adj[v]:
                if w in Rs:
                    if idx[w] > idx[v]:
                        for c in range(4):
                            s.add_clause([-var(v, c), -var(w, c)])
                else:
                    if col[w] >= 0:
                        s.add_clause([-var(v, col[w])])
        # phases from old colours
        s.set_phases([var(v, col[v]) for v in R if col[v] >= 0])
        t = time.time()
        res = solve_limited(s, tlimit=tl)
        dt = time.time() - t
        log(f"    lns r={r}: free {len(R)} vertices -> {res} in {dt:.1f}s")
        if res:
            m = set(l for l in s.get_model() if l > 0)
            out = list(col)
            for v in R:
                out[v] = next(c for c in range(4) if var(v, c) in m)
            s.delete()
            # outside vertices not in R must all be coloured
            if any(c < 0 for c in out):
                return None
            return out
        s.delete()
        if r >= rmax or res is None:
            return None
        # grow region by one BFS layer
        r += 1
        nf = []
        for v in frontier:
            for w in adj[v]:
                if w not in dist:
                    dist[w] = r
                    nf.append(w)
        frontier = nf
        if not nf and len(dist) >= n:
            return None
