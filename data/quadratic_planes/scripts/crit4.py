"""crit4.py p k IN.npy OUT [K]: shrink a non-K-colourable induced subgraph of the level-k p-adic plane
G_k = Cay((Z/p^k)^2, T_k) (vertex codes x * p^k + y in IN.npy) to a vertex-critical one, with one incremental
CaDiCaL solver: a selector per vertex switches on its 'some colour' clause; for each vertex v (lowest degree
first) solve without v: SAT means v is necessary, UNSAT gives a core of selectors and the set shrinks to it.
Saves OUT.npy (codes) and prints the size. No symmetry breaking (so the cores are honest subsets)."""
import sys, time
import numpy as np
from collections import defaultdict
from pysat.solvers import Solver
p, k, IN, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
K = int(sys.argv[5]) if len(sys.argv) > 5 else 4
m = p ** k; t0 = time.time()
roots = defaultdict(list)
for y in range(m):
    roots[(y * y) % m].append(y)
T = np.array([(a, b) for a in range(m) for b in roots.get((1 - a * a) % m, [])], dtype=np.int64)
V = np.unique(np.load(IN).astype(np.int64)); n = len(V)
pos = -np.ones(m * m, dtype=np.int64); pos[V] = np.arange(n)
vx, vy = V // m, V % m
adj = [set() for _ in range(n)]
for i in range(0, len(T), 64):
    tb = T[i:i + 64]
    w = pos[(((vx[:, None] + tb[None, :, 0]) % m) * m + (vy[:, None] + tb[None, :, 1]) % m)]
    for a, row in enumerate(w.tolist()):
        for b in row:
            if b >= 0:
                adj[a].add(b)
E = [(a, b) for a in range(n) for b in adj[a] if a < b]
print(f"p={p} level {k}: {n} vertices, {len(E)} edges  [{time.time() - t0:.1f}s]", flush=True)
x = lambda v, c: K * v + c + 1
sel = lambda v: K * n + v + 1
s = Solver(name="cadical153")
for v in range(n):
    s.add_clause([-sel(v)] + [x(v, c) for c in range(K)])
for a, b in E:
    for c in range(K):
        s.add_clause([-x(a, c), -x(b, c)])
PIN = []
for a, b in E:
    c3 = adj[a] & adj[b]
    if c3:
        PIN = [a, b, min(c3)]; break
for c, v in enumerate(PIN):                 # colours of a triangle fixed: no loss for K >= 3
    s.add_clause([x(v, c)])
print(f"pinned triangle {PIN}", flush=True)


def kcore(S):
    S = set(S)
    stack = [v for v in S if len(adj[v] & S) < K and v not in PIN]
    while stack:
        v = stack.pop()
        if v not in S or v in PIN:
            continue
        S.discard(v)
        for w in adj[v]:
            if w in S and len(adj[w] & S) < K:
                stack.append(w)
    return S


def refute(S):
    S = set(S) | set(PIN)
    if s.solve(assumptions=[sel(v) for v in sorted(S)]):
        return None
    return {l - K * n - 1 for l in s.get_core() if l > K * n} | set(PIN)


S = kcore(range(n))
C = refute(S)
assert C is not None, "the input is K-colourable"
S = kcore(C)
print(f"first core: {len(S)}  [{time.time() - t0:.0f}s]", flush=True)
necessary = set(PIN); tests = 0
for v in sorted(S, key=lambda v: (len(adj[v] & S), v)):
    if v not in S or v in necessary:
        continue
    tests += 1
    C = refute(kcore(S - {v}))
    if C is None:
        necessary.add(v)
    else:
        S = kcore(C)
    if tests % 50 == 0:
        print(f"  {tests} tests: {len(S)} vertices, {len(necessary)} necessary  [{time.time() - t0:.0f}s]", flush=True)
assert necessary >= S
np.save(out + ".npy", V[sorted(S)])
print(f"vertex-critical: {len(S)} vertices, {sum(len(adj[v] & S) for v in S) // 2} edges; saved {out}.npy  [{time.time() - t0:.0f}s]", flush=True)
