"""Deletion-based minimisation of a chi_c = p/q witness (lazy SAT with one selector per vertex).

usage: python3 minimise.py WITNESS.json OUT.json [conflict_budget] [order]
WITNESS.json: {"d", "denominator", "points", "edges", "colouring", "cycles", ...} (as written by grow.py); order:
"degree" (lowest degree first) or an integer seed for a random order.  The value p/q is read from the file
(keys "p", "q"; 7/2 if absent); below, read (p,q), q mod p and x_v0..x_v(p-1) for general p/q.
Vertex v is kept iff its selector a_v is assumed; "v has a colour" is the clause (-a_v, x_v0..x_v6), edge, tightness
and cycle clauses are unconditional (an uncoloured vertex satisfies them).  For each v in turn: solve with the kept
vertices minus v; SAT -> list tight cycles of the colouring inside the kept set minus v and repeat; acyclic -> v is
needed; UNSAT -> drop v and every kept vertex outside the core.  A call over the conflict budget keeps v (the
first and the last call, on the whole kept set, have no budget).
The colour of one kept vertex (the root) is fixed to 0 under its own selector, as rotating colours keeps tightness."""
import json, sys, time, random
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
BUD = int(sys.argv[3]) if len(sys.argv) > 3 else 3000000
ORDER = sys.argv[4] if len(sys.argv) > 4 else 'degree'
PP, QQ = W.get('p', 7), W.get('q', 2)
FD = sorted({dl % PP for dl in range(1 - QQ, QQ)})    # forbidden colour differences on an edge
n = len(W['points']); E = [tuple(e) for e in W['edges']]
adj = [[] for _ in range(n)]
for i, j in E:
    adj[i].append(j); adj[j].append(i)
x = lambda v, c: PP * v + c + 1
a = lambda v: PP * n + v + 1
f = lambda v: (PP + 1) * n + v + 1   # root selectors
nv = [(PP + 2) * n]
s = Solver(name='cadical153')
for v in range(n):
    s.add_clause([-a(v)] + [x(v, c) for c in range(PP)])
    s.add_clause([-f(v), x(v, 0)])
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
cycles = [list(c) for c in W['cycles']]
for c in cycles:
    s.add_clause([-T(c[k], c[(k + 1) % len(c)]) for k in range(len(c))])
kept = set(range(n))
t0 = time.time()
def test(drop, budget=True):
    """True if the kept set minus `drop` still has no admissible colouring (UNSAT); also returns the core."""
    K = kept - set(drop)
    root = min(K)
    sub_adj = None
    while True:
        assum = [a(v) for v in sorted(K)] + [f(root)]
        if budget:
            s.conf_budget(BUD)
            r = s.solve_limited(assumptions=assum)
        else:
            r = s.solve(assumptions=assum)
        if r is None:
            return None, None
        if not r:
            cl = s.get_core() or []
            core = set(l - PP * n - 1 for l in cl if PP * n < l <= (PP + 1) * n)
            if f(root) in cl:
                core.add(root)          # the fixed colour of the root is part of the refutation
            return True, core
        m = s.get_model(); ms = set(l for l in m if l > 0)
        col = [next((c for c in range(PP) if x(v, c) in ms), -1) if v in K else -1 for v in range(n)]
        if sub_adj is None:
            sub_adj = [[w for w in adj[v] if w in K] if v in K else [] for v in range(n)]
        # vertices outside K get isolated; colour -1 vertices have no arcs
        colx = [c if c >= 0 else 0 for c in col]
        cyc = tight_cycles(n, sub_adj, colx, 40)
        cyc = [c for c in cyc if all(v in K for v in c)]
        if not cyc:
            return False, None
        for c in cyc:
            s.add_clause([-T(c[k], c[(k + 1) % len(c)]) for k in range(len(c))])
            cycles.append(c)
deg = [len(adj[v]) for v in range(n)]
order = sorted(range(n), key=lambda v: (deg[v], v)) if ORDER == 'degree' else random.Random(int(ORDER)).sample(range(n), n)
ok, core = test([], budget=False)
assert ok, 'the witness is not UNSAT'
print(f'start: {n} vertices; core {len(core)} [{time.time() - t0:.0f}s]', flush=True)
kept &= core
k = 0
for v in order:
    if v not in kept:
        continue
    k += 1
    ok, core = test([v])
    if ok:
        kept.discard(v); kept &= core
    if k % 10 == 0 or ok is None:
        print(f'{k}: v={v} -> {"drop" if ok else ("budget" if ok is None else "keep")}; kept {len(kept)} [{time.time() - t0:.0f}s]', flush=True)
    if k % 25 == 0:
        json.dump(sorted(kept), open(outp + '.partial', 'w'))
ok, core = test([], budget=False)
assert ok
print(f'final: {len(kept)} vertices [{time.time() - t0:.0f}s]', flush=True)
K = sorted(kept); ren = {v: i for i, v in enumerate(K)}
E2 = [[ren[i], ren[j]] for i, j in E if i in kept and j in kept]
cyc2 = [[ren[v] for v in c] for c in cycles if all(v in kept for v in c)]
pq = {} if (PP, QQ) == (7, 2) else {'p': PP, 'q': QQ}
json.dump({'d': W['d'], 'denominator': W['denominator'], **pq, 'seed': 'minimise(' + sys.argv[1] + ')',
           'points': [W['points'][v] for v in K], 'edges': E2, 'colouring': [W['colouring'][v] for v in K],
           'cycles': cyc2}, open(outp, 'w'))
print('wrote', outp, len(K), 'points', len(E2), 'edges', len(cyc2), 'cycles', flush=True)
