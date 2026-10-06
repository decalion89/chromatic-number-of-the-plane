"""per_build.py CORE.json POOL.json [POOL2.json ...] OUT.json: a finite point set H in which the base-point periods of
the core relations are tied to the 23 basis relations by homotopies (Lemma P), built as follows.

Walk of a relation vector v (27 integers): from 0, |v_u| steps of sign(v_u) u for u = 0..26 in ascending order.  A
chain from 0 to each core relation adds or subtracts one basis relation r_j at a time, staying inside the pool of
short relations (BFS).  Each chain step rho' = rho + s r_j is a triple: the walk W_rho W_{r_j}^s is bubble-sorted into
W_{rho'} (a swap of adjacent steps adds the fourth corner of a square to H; adjacent opposite steps cancel).  So
per(W_{rho'}) = per(W_rho) + s per(W_{r_j}) for every 4-colouring of H without tight cycles.  H = all vertices met.
Then CP-SAT checks that the integer system (one variable per relation met; the open range of Lemma P for each; one
equation per triple) is infeasible, which proves chi_c(H) >= 4 if Lemma P holds; the SAT check of H is separate."""
import json, sys, time, gzip, os
from fractions import Fraction as Fr
from math import floor, ceil
from collections import deque
from ortools.sat.python import cp_model

t0 = time.time()
C = json.load(gzip.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "at_four", "cert311_open_4.json.gz"), "rt")); U = C["units"]; Rb = [tuple(r) for r in C["relations"]]
m = len(U); dim = len(U[0])
core = [tuple(r["rho"]) for r in json.load(open(sys.argv[1]))]
pool = set()
for f in sys.argv[2:-1]:
    for r in json.load(open(f)):
        pool.add(tuple(r["rho"])); pool.add(tuple(-x for x in r["rho"]))
for r in Rb:
    pool.add(r); pool.add(tuple(-x for x in r))
out = sys.argv[-1]
Z = (0,) * m
pool.add(Z)
add = lambda a, b, s: tuple(x + s * y for x, y in zip(a, b))
# BFS from 0 over the pool with moves +- r_j
par = {Z: None}; q = deque([Z])
while q:
    v = q.popleft()
    for j, r in enumerate(Rb):
        for s in (1, -1):
            w = add(v, r, s)
            if w in pool and w not in par:
                par[w] = (v, j, s); q.append(w)
miss = [c for c in core if c not in par]
print(f"pool {len(pool)}, reached {len(par)}; core {len(core)}, unreached {len(miss)} [{time.time()-t0:.1f}s]",
      flush=True)
# unreached: best-first search from the relation to the reached set, intermediate relations of l1 <= LMAX
import heapq
LMAX = int(__import__("os").environ.get("LMAX", "22"))
l1 = lambda v: sum(abs(a) for a in v)
for c in miss:
    if c in par:
        continue
    prev = {c: None}; h = [(l1(c), c)]; hit = None; n_exp = 0
    while h and n_exp < 400000:
        _, v = heapq.heappop(h); n_exp += 1
        if v in par:
            hit = v; break
        for j, r in enumerate(Rb):
            for s in (1, -1):
                w = add(v, r, s)                       # v = w - s r_j, i.e. w = v + s r_j
                if w in prev or l1(w) > LMAX:
                    continue
                prev[w] = (v, j, s)
                heapq.heappush(h, (l1(w), w))
    assert hit is not None, ("no chain for", c)
    # hit is reached; walk back from hit to c: each step v -> w = v + s r_j means v = w + (-s) r_j
    w = hit
    while prev[w] is not None:
        v, j, s = prev[w]
        par[v] = (w, j, -s)
        w = v
    print(f"  chain for a relation of l1 {l1(c)}: {n_exp} expansions", flush=True)
# triples along the chains
rels = {}; triples = []
def rid(v):
    if v not in rels:
        rels[v] = len(rels)
    return rels[v]
rid(Z)
for j, r in enumerate(Rb):
    rid(r)
done = set()
for c in core:
    v = c
    while par[v] is not None:
        if v in done:
            break
        done.add(v)
        u, j, s = par[v]
        triples.append((rid(u), s, rid(Rb[j]), rid(v)))       # W_u W_{r_j}^s ~ W_v
        v = u
print(f"{len(rels)} relations, {len(triples)} triples", flush=True)
R = [None] * len(rels)
for v, i in rels.items():
    R[i] = v


def steps_of(v, sign=1):
    st = []
    for j, c in enumerate(v):
        st += [(j, 1 if c > 0 else -1)] * abs(c)
    if sign < 0:
        st = [(j, -s) for (j, s) in reversed(st)]
    return st


def points_of(steps):
    p = [0] * dim; pts = [tuple(p)]
    for j, s in steps:
        for t in range(dim):
            p[t] += s * U[j][t]
        pts.append(tuple(p))
    return pts


H = set(); swaps = canc = 0
for v in R:
    H.update(points_of(steps_of(v)))
for (i, s, j, k) in triples:
    st = steps_of(R[i]) + steps_of(R[j], s)
    pts = points_of(st); H.update(pts)
    changed = True
    while changed:
        changed = False
        a = 0
        while a < len(st) - 1:
            x, y = st[a], st[a + 1]
            if x[0] == y[0] and x[1] == -y[1]:
                del st[a:a + 2]; del pts[a + 1:a + 3]; canc += 1; changed = True; a = max(a - 1, 0); continue
            if x[0] > y[0]:
                base = pts[a]
                new = tuple(base[t] + y[1] * U[y[0]][t] for t in range(dim))
                H.add(new); st[a], st[a + 1] = y, x; pts[a + 1] = new; swaps += 1; changed = True
            a += 1
    assert st == steps_of(R[k]), "chain step does not sort to its target"
print(f"|H| = {len(H)} points, {swaps} swaps, {canc} cancellations [{time.time()-t0:.1f}s]", flush=True)
# the integer system on the relations met
md = cp_model.CpModel()
x = [md.NewIntVar(-500, 500, f"x{i}") for i in range(len(R))]
md.Add(x[0] == 0)
for i, v in enumerate(R):
    if i == 0:
        continue
    P = sum(a for a in v if a > 0); N = -sum(a for a in v if a < 0)
    md.Add(x[i] >= floor(Fr(P - 3 * N, 4)) + 1); md.Add(x[i] <= ceil(Fr(3 * P - N, 4)) - 1)
for (i, s, j, k) in triples:
    md.Add(x[k] == x[i] + s * x[j])
sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = 600; sv.parameters.num_workers = 1
st = sv.Solve(md)
print("period system on H:", sv.StatusName(st), f"[{time.time()-t0:.1f}s]", flush=True)
json.dump({"D": C["D"], "units": U, "points": sorted(H), "relations": [list(v) for v in R],
           "triples": triples}, open(out, "w"))
