"""Colouring-guided growth of a finite unit-distance graph H over Q(sqrt d) with chi_c(H) = p/q (by default 7/2).

usage: python3 grow.py GRAPH.json SEED OUT_PREFIX [max_points] [add_per_round] [p q]
  GRAPH.json: {"d", "D", "points"} (data/quadratic_planes/q<d>.json); SEED: "A" (its points) or "A+A" (their sumset).
  p q: the value p/q, 2 <= 2q <= p (default 7 2).  Below, for general (p, q), read (p,q) for (7,2), q mod p for
  2 mod 7, and the differences in (-q, q) mod p for {0, 1, -1} mod 7; the output then also records p and q.
Lazy SAT: look for a (7,2)-colouring of H in which every listed directed cycle has a non-tight arc (an arc a -> b is
tight when c(b) - c(a) = 2 mod 7).  When the colouring found has tight cycles, they are listed and the search goes
on; when its tight digraph is acyclic, H grows: the candidates p = x + u (x in H, u a unit vector with denominator D)
whose unit neighbours in H leave p no colour (every colour k has a neighbour y with k - c(y) in {0, 1, -1} mod 7) are
added, most neighbours first, at most add_per_round of them; if there is none, the candidates with one colour left.
Stops when no such colouring exists: then every (7,2)-colouring of H has a tight cycle, so chi_c(H) >= 7/2 by
Guichard's lemma (Lemma 20 of papers/three-colours), and the (7,2)-colouring written out shows chi_c(H) = 7/2; or at
max_points.  Writes OUT_PREFIX.json: d, denominator, points, edges (all unit pairs), a (7,2)-colouring, cycles.
The colour of vertex 0 is fixed to 0 (rotating the colours keeps tightness).  With python-sat 1.9 (CaDiCaL 1.5.3),
python3 grow.py data/quadratic_planes/q11.json A OUT 3000 200 gives the 653-vertex graph of witness_q11grow, and
python3 grow.py q7_seed.json A OUT 2000 200 3 1 the 607-vertex graph from which witness_q7 was cut."""
import json, math, sys, time
import numpy as np
from pysat.solvers import Solver


def unit_vectors(d, D):
    out = set()
    B = math.isqrt(D * D // d)
    for b in range(-B, B + 1):
        for e in range(-B, B + 1):
            M = D * D - d * (b * b + e * e)
            if M < 0:
                continue
            if b == 0 and e == 0:
                for a in range(-D, D + 1):
                    c2 = M - a * a
                    if c2 >= 0:
                        c = math.isqrt(c2)
                        if c * c == c2:
                            out.add((a, 0, c, 0)); out.add((a, 0, -c, 0))
                continue
            g = math.gcd(b, e); bp, ep = b // g, e // g
            nn = bp * bp + ep * ep
            if M % nn:
                continue
            k = math.isqrt(M // nn)
            if k * k * nn != M:
                continue
            for kk in {k, -k}:
                out.add((-kk * ep, b, kk * bp, e))
    for (a, b, c, e) in out:
        assert a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0
    return sorted(out)


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


path, seed, outp = sys.argv[1], sys.argv[2], sys.argv[3]
MAXP = int(sys.argv[4]) if len(sys.argv) > 4 else 60000
ADD = int(sys.argv[5]) if len(sys.argv) > 5 else 400
PP, QQ = (int(sys.argv[6]), int(sys.argv[7])) if len(sys.argv) > 7 else (7, 2)
assert 2 <= 2 * QQ <= PP
FD = sorted({dl % PP for dl in range(1 - QQ, QQ)})    # forbidden colour differences on an edge: [0, 1, 6] for (7,2)
g = json.load(open(path)); d, D = g['d'], g['D']
A = [tuple(p) for p in g['points']]
P0 = sorted(set(A)) if seed == 'A' else sorted({tuple(x + y for x, y in zip(p, q)) for p in A for q in A})
U = unit_vectors(d, D)
OFF, M = 1 << 20, 1 << 21          # coordinates must stay in (-OFF, OFF)
def key(p): return ((((p[0] + OFF) * M + p[1] + OFF) * M + p[2] + OFF) * M + p[3] + OFF)
def unkey(k):
    e = k % M - OFF; k //= M; c = k % M - OFF; k //= M; b = k % M - OFF; k //= M; a = k - OFF
    return (int(a), int(b), int(c), int(e))
UD = np.array([(((u[0]) * M + u[1]) * M + u[2]) * M + u[3] for u in U], dtype=object)
t0 = time.time()
pts = []; idx = {}; adj = []
s = Solver(name='cadical153'); nv = [0]
base = []
def newvar():
    nv[0] += 1; return nv[0]
def add_point(p):
    v = len(pts); pts.append(p); k = key(p); idx[k] = v; adj.append([])
    b = nv[0] + 1; nv[0] += PP; base.append(b)
    s.add_clause([b + c for c in range(PP)])
    for du in UD:
        w = idx.get(k + du)
        if w is not None and w != v:
            adj[v].append(w); adj[w].append(v)
            bw = base[w]
            for c in range(PP):
                for dl in FD:
                    s.add_clause([-(b + c), -(bw + (c + dl) % PP)])
for p in P0:
    add_point(p)
if len(pts) > 0:
    s.add_clause([base[0]])          # colour of vertex 0 fixed to 0 (rotating colours preserves tightness)
tv = {}; cycles = []
def T(a, b):
    if (a, b) not in tv:
        t = newvar(); tv[(a, b)] = t
        for c in range(PP):
            s.add_clause([-(base[a] + c), -(base[b] + (c + QQ) % PP), t])
    return tv[(a, b)]
F = [sum(1 << ((c + dl) % PP) for dl in FD) for c in range(PP)]
FULL = (1 << PP) - 1
ne = lambda: sum(len(a) for a in adj) // 2
print(f'd={d} D={D} seed={seed}: {len(pts)} points, {ne()} unit pairs, {len(U)} unit vectors [{time.time() - t0:.0f}s]', flush=True)
rnd = 0
while True:
    rnd += 1
    if not s.solve():
        print(f'round {rnd}: UNSAT with {len(pts)} points, {ne()} unit pairs, {len(cycles)} cycles: chi_c(H) = {PP}/{QQ} [{time.time() - t0:.0f}s]', flush=True)
        s2 = Solver(name='cadical153')
        for v in range(len(pts)):
            s2.add_clause([base[v] + c for c in range(PP)])
            for w in adj[v]:
                if w > v:
                    for c in range(PP):
                        for dl in FD:
                            s2.add_clause([-(base[v] + c), -(base[w] + (c + dl) % PP)])
        assert s2.solve(), f'no ({PP},{QQ})-colouring at all'
        ms = set(l for l in s2.get_model() if l > 0)
        col = [next(c for c in range(PP) if base[v] + c in ms) for v in range(len(pts))]
        E = sorted((v, w) for v in range(len(pts)) for w in adj[v] if w > v)
        pq = {} if (PP, QQ) == (7, 2) else {'p': PP, 'q': QQ}
        json.dump({'d': d, 'denominator': D, **pq, 'seed': seed, 'source': path, 'points': [list(p) for p in pts],
                   'edges': [list(e) for e in E], 'colouring': col, 'cycles': cycles}, open(outp + '.json', 'w'))
        print('wrote', outp + '.json', flush=True)
        break
    ms = set(l for l in s.get_model() if l > 0)
    col = [next(c for c in range(PP) if base[v] + c in ms) for v in range(len(pts))]
    cyc = tight_cycles(len(pts), adj, col, 40)
    if cyc:
        for c in cyc:
            s.add_clause([-T(c[k], c[(k + 1) % len(c)]) for k in range(len(c))])
            cycles.append(c)
        if rnd % 10 == 0:
            print(f'round {rnd}: +{len(cyc)} cycles (total {len(cycles)}), {len(pts)} points [{time.time() - t0:.0f}s]', flush=True)
        continue
    if len(pts) >= MAXP:
        print(f'round {rnd}: acyclic tight digraph and {len(pts)} points >= max; stop [{time.time() - t0:.0f}s]', flush=True)
        break
    # growth: candidates p = x + u
    hk = np.array([key(p) for p in pts], dtype=object)
    cand = {}
    for x, kx in enumerate(hk):
        for du in UD:
            kp = kx + du
            if kp in idx:
                continue
            cand.setdefault(kp, 0)
    ck = list(cand.keys())
    forb = [0] * len(ck); nbrs = [0] * len(ck)
    for i, kp in enumerate(ck):
        f = 0; nb = 0
        for du in UD:
            w = idx.get(kp + du)
            if w is not None:
                f |= F[col[w]]; nb += 1
        forb[i] = f; nbrs[i] = nb
    blocked = [i for i in range(len(ck)) if forb[i] == FULL]
    if blocked:
        blocked.sort(key=lambda i: -nbrs[i]); pick = blocked[:ADD]; kind = 'blocked'
    else:
        one = [i for i in range(len(ck)) if bin(FULL & ~forb[i]).count('1') == 1]
        one.sort(key=lambda i: -nbrs[i]); pick = one[:ADD]; kind = 'one colour left'
    if not pick:
        print(f'round {rnd}: acyclic tight digraph and no candidate; stop [{time.time() - t0:.0f}s]', flush=True)
        break
    for i in pick:
        add_point(unkey(ck[i]))
    print(f'round {rnd}: grew by {len(pick)} ({kind}; {len(blocked)} blocked of {len(ck)} candidates) -> {len(pts)} points, {ne()} unit pairs [{time.time() - t0:.0f}s]', flush=True)
