"""Vertex-criticality certificates for a chi_c = p/q witness (p/q read from the file, 7/2 by default): for every
vertex v, a (p,q)-colouring of H - v whose tight digraph is acyclic (so chi_c(H - v) < p/q).  If some H - v still has no such colouring, v is reported (the
witness is then not vertex-critical; minimise again).

usage: python3 critical.py WITNESS.json OUT.json   (check the output with check_critical.py, for 7/2)"""
import json, sys, time
from pysat.solvers import Solver


def tight_cycles(n, adj, col, limit):
    """Directed cycles of the tight digraph (arcs a->b with col[b]-col[a] = q mod p): Tarjan SCCs, then from one
    vertex of each nontrivial SCC a shortest cycle through it (BFS inside the SCC); at most `limit` cycles."""
    out = [[b for b in adj[a] if (col[b] - col[a]) % PP == QQ] for a in range(n)]
    index = [None] * n; low = [0] * n; onst = [False] * n; st = []; comp = [-1] * n; ncomp = 0; t = 0
    for s in range(n):
        if index[s] is not None:
            continue
        work = [(s, 0)]
        while work:
            v, k = work.pop()
            if k == 0:
                index[v] = low[v] = t; t += 1; st.append(v); onst[v] = True
            recurse = False
            for kk in range(k, len(out[v])):
                w = out[v][kk]
                if index[w] is None:
                    work.append((v, kk + 1)); work.append((w, 0)); recurse = True; break
                elif onst[w]:
                    low[v] = min(low[v], index[w])
            if recurse:
                continue
            if low[v] == index[v]:
                while True:
                    w = st.pop(); onst[w] = False; comp[w] = ncomp
                    if w == v:
                        break
                ncomp += 1
            if work:
                u = work[-1][0]
                low[u] = min(low[u], low[v])
    members = {}
    for v in range(n):
        members.setdefault(comp[v], []).append(v)
    cycles = []
    for cid, vs in members.items():
        if len(vs) < 2:
            continue
        for s in vs[:max(1, limit // max(1, len(members)))]:
            if len(cycles) >= limit:
                break
            prev = {s: None}; frontier = [s]; found = None
            while frontier and found is None:
                nxt = []
                for v in frontier:
                    for w in out[v]:
                        if comp[w] != cid:
                            continue
                        if w == s:
                            found = v; break
                        if w not in prev:
                            prev[w] = v; nxt.append(w)
                    if found is not None:
                        break
                frontier = nxt
            if found is None:
                continue
            cyc = []; v = found
            while v is not None:
                cyc.append(v); v = prev[v]
            cycles.append(cyc[::-1])
        if len(cycles) >= limit:
            break
    return cycles


W = json.load(open(sys.argv[1])); outp = sys.argv[2]
PP, QQ = W.get('p', 7), W.get('q', 2)
FD = sorted({dl % PP for dl in range(1 - QQ, QQ)})    # forbidden colour differences on an edge
n = len(W['points']); E = [tuple(e) for e in W['edges']]
adj = [[] for _ in range(n)]
for i, j in E:
    adj[i].append(j); adj[j].append(i)
x = lambda v, c: PP * v + c + 1
a = lambda v: PP * n + v + 1
nv = [(PP + 1) * n]
s = Solver(name='cadical153')
for v in range(n):
    s.add_clause([-a(v)] + [x(v, c) for c in range(PP)])
for i, j in E:
    for c in range(PP):
        for dl in FD:
            s.add_clause([-x(i, c), -x(j, (c + dl) % PP)])
tv = {}
def T(p, q):
    if (p, q) not in tv:
        nv[0] += 1; tv[(p, q)] = nv[0]
        for c in range(PP):
            s.add_clause([-x(p, c), -x(q, (c + QQ) % PP), nv[0]])
    return tv[(p, q)]
for c in W['cycles']:
    s.add_clause([-T(c[k], c[(k + 1) % len(c)]) for k in range(len(c))])
t0 = time.time(); certs = []; bad = []
for v in range(n):
    K = set(range(n)) - {v}
    sub = [[w for w in adj[u] if w != v] if u != v else [] for u in range(n)]
    while True:
        if not s.solve(assumptions=[a(u) for u in sorted(K)]):
            bad.append(v); break
        ms = set(l for l in s.get_model() if l > 0)
        col = [next((c for c in range(PP) if x(u, c) in ms), 0) if u != v else 0 for u in range(n)]
        cyc = [c for c in tight_cycles(n, sub, col, 40) if v not in c]
        if not cyc:
            col[v] = -1
            certs.append(col); break
        for c in cyc:
            s.add_clause([-T(c[k], c[(k + 1) % len(c)]) for k in range(len(c))])
    if v % 25 == 0:
        print(f'{v}/{n} [{time.time() - t0:.0f}s]', flush=True)
print(f'done: {len(certs)} certificates, not critical at {bad} [{time.time() - t0:.0f}s]', flush=True)
json.dump({'critical_colourings': certs, 'not_critical': bad}, open(outp, 'w'))
