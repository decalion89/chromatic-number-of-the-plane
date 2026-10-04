#!/usr/bin/env python3
"""Cycle-list-independent check of claim 3 (referee code).

Usage:  python3 referee_cegar.py W.json.gz OUT.cnf OUT_cycles.json [FIXED_VERTEX] [SOLVER]

Ignores the cycles stored in the witness.  Starting from an EMPTY list, repeatedly asks a SAT solver (pysat) for a
(7,2)-colouring of H in which every cycle found so far has a non-tight arc (colour of FIXED_VERTEX fixed to 0;
default: vertex 0, not the repository's fixed vertex).  If the answer's tight digraph is acyclic, that colouring
would be a COUNTEREXAMPLE (chi_c(H) < 7/2) and is reported.  Otherwise, for every nontrivial strongly connected
component of the tight digraph, a shortest directed cycle through some vertex of it is added to the list.  When the
solver says UNSAT, the final formula (E1 layout, own cycle list) is written to OUT.cnf for kissat/drat-trim.
"""
import gzip, json, sys, time, hashlib
from collections import defaultdict, deque
from pysat.solvers import Solver


def tarjan(nodes, out):
    index, low, onst, st, comps = {}, {}, set(), [], []
    counter = [0]
    sys.setrecursionlimit(100000)

    def strong(v):
        index[v] = low[v] = counter[0]; counter[0] += 1
        st.append(v); onst.add(v)
        for w in out[v]:
            if w not in index:
                strong(w); low[v] = min(low[v], low[w])
            elif w in onst:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = st.pop(); onst.discard(w); comp.append(w)
                if w == v:
                    break
            comps.append(comp)
    for v in nodes:
        if v not in index:
            strong(v)
    return comps


def shortest_cycle_through(s, out, allowed):
    # BFS from s inside 'allowed' back to s
    prev = {s: None}
    dq = deque([s])
    while dq:
        u = dq.popleft()
        for w in out[u]:
            if w not in allowed:
                continue
            if w == s:
                path = [u]
                while prev[path[-1]] is not None:
                    path.append(prev[path[-1]])
                return path[::-1]
            if w not in prev:
                prev[w] = u
                dq.append(w)
    return None


if not __debug__:
    sys.exit('run without -O: this script relies on assert statements')


def main():
    with gzip.open(sys.argv[1], 'rt') as f:
        W = json.load(f)
    out_cnf, out_cyc = sys.argv[2], sys.argv[3]
    fv = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    sname = sys.argv[5] if len(sys.argv) > 5 else 'cadical153'
    n = len(W['points'])
    E = [tuple(e) for e in W['edges']]
    x = lambda v, k: k * n + v + 1
    base = []
    for v in range(n):
        base.append([x(v, k) for k in range(7)])
        base += [[-x(v, k), -x(v, l)] for k in range(7) for l in range(k + 1, 7)]
    for (u, w) in E:
        for k in range(7):
            for dl in (0, 1, 6):
                base.append([-x(u, k), -x(w, (k + dl) % 7)])
    base.append([x(fv, 0)])
    arcvar = {}
    nxt = [7 * n]
    extra = []

    def N(a):
        if a not in arcvar:
            nxt[0] += 1
            arcvar[a] = nxt[0]
            for k in range(7):
                extra.append([-nxt[0], -x(a[0], k), -x(a[1], (k + 2) % 7)])
        return arcvar[a]

    s = Solver(name=sname, bootstrap_with=base)
    cycles = []
    seen = set()
    t0 = time.time()
    it = 0
    while True:
        it += 1
        if not s.solve():
            break
        model = s.get_model()
        tv = set(l for l in model if l > 0)
        col = [None] * n
        for v in range(n):
            ks = [k for k in range(7) if x(v, k) in tv]
            assert len(ks) == 1
            col[v] = ks[0]
        for (u, w) in E:
            assert (col[w] - col[u]) % 7 in (2, 3, 4, 5)
        out = defaultdict(list)
        for (u, w) in E:
            dd = (col[w] - col[u]) % 7
            if dd == 2:
                out[u].append(w)
            elif dd == 5:
                out[w].append(u)
        comps = [c for c in tarjan(range(n), out) if len(c) > 1]
        if not comps:
            print(f"COUNTEREXAMPLE: (7,2)-colouring with acyclic tight digraph found at iteration {it}: {col}")
            sys.exit(2)
        added = 0
        for comp in comps:
            cs = set(comp)
            # several cycles per component: through up to 3 vertices of it
            for sv in comp[:3]:
                C = shortest_cycle_through(sv, out, cs)
                key = min(tuple(C[i:] + C[:i]) for i in range(len(C)))
                if key in seen:
                    continue
                seen.add(key)
                cycles.append(list(C))
                m = len(C)
                lits = [N((C[i], C[(i + 1) % m])) for i in range(m)]
                for c in extra:
                    s.add_clause(c)
                extra.clear()
                s.add_clause(lits)
                added += 1
        if it % 50 == 0:
            print(f"  iteration {it}: {len(cycles)} cycles, {time.time()-t0:.1f} s", flush=True)
    el = time.time() - t0
    print(f"UNSAT after {it} iterations, {len(cycles)} cycles of lengths "
          f"{sorted(set(len(C) for C in cycles))}, {el:.1f} s (solver {sname}, fixed vertex {fv})")
    # final formula, E1 layout with own cycles
    cls = list(base)
    arcl = {}
    nv = 7 * n
    for C in cycles:
        m = len(C)
        for i in range(m):
            a = (C[i], C[(i + 1) % m])
            if a not in arcl:
                nv += 1
                arcl[a] = nv
                for k in range(7):
                    cls.append([-nv, -x(a[0], k), -x(a[1], (k + 2) % 7)])
    for C in cycles:
        m = len(C)
        cls.append([arcl[(C[i], C[(i + 1) % m])] for i in range(m)])
    txt = f"p cnf {nv} {len(cls)}\n" + ''.join(' '.join(map(str, c)) + ' 0\n' for c in cls)
    open(out_cnf, 'w').write(txt)
    json.dump(cycles, open(out_cyc, 'w'))
    print(f"wrote {out_cnf}: {nv} variables, {len(cls)} clauses, sha256 {hashlib.sha256(txt.encode()).hexdigest()}")


if __name__ == '__main__':
    main()
