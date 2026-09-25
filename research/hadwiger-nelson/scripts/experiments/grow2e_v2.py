"""Gadget growth, v2: symmetry broken, checkpointed, resumable.

Same loop as grow2e -- ask for a 5-colouring with c(a) = c(b), insert the pool
points whose neighbourhoods that colouring shows all five colours on -- but
(1) a triangle of the graph is pinned to colours 0, 1, 2, which any colouring
can be permuted into, so the solver never re-explores a colour relabelling;
(2) the graph is written every 5 iterations so a stopped run resumes where it
was; (3) selection mixes richness with nearness to the pair.

usage: grow2e_v2.py <in.json with points, A, B> <out.json> [R] [near_weight]
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
t0 = time.time(); K = 5
IN, OUT = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 20
WN = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
d = json.load(open(IN))
F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
A, Bi = d["A"], d["B"]
assert V[A].dist2(V[Bi]) == F.rational(4)
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
vkey = {key(p): i for i, p in enumerate(V)}
G = build_graph(V)
U = {}
for a, b in G.edges():
    for u in (V[b] - V[a], V[a] - V[b]): U.setdefault(key(u), u)
U = list(U.values())
pool = {}
def feed(i):
    for u in U:
        p = V[i] + u; k = key(p)
        if k in vkey: continue
        e = pool.get(k)
        if e is None: pool[k] = [p, {i}]
        else: e[1].add(i)
for i in range(len(V)): feed(i)
def save(tag):
    json.dump({"field_generators": list(F.gens), "A": A, "B": Bi, "status": tag,
               "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in V]},
              open(OUT, "w"))
print(f"{IN}: n={G.n} m={G.m}, {len(U)} directions, pool {len(pool)}; pair {A},{Bi}   [{time.time()-t0:.0f}s]", flush=True)
ax, ay, bx, by = V[A].fx, V[A].fy, V[Bi].fx, V[Bi].fy
near = lambda p: min(((p.fx - ax) ** 2 + (p.fy - ay) ** 2) ** 0.5, ((p.fx - bx) ** 2 + (p.fy - by) ** 2) ** 0.5)
for it in range(1, 100000):
    n = G.n; E = list(G.edges())
    adj = [set() for _ in range(n)]
    for a, b in E: adj[a].add(b); adj[b].add(a)
    tri = next(((a, b, c) for a, b in E for c in adj[a] & adj[b]), None)
    X = lambda v, c: 1 + v * K + c
    s = Solver(name="cd19")
    for v in range(n): s.add_clause([X(v, c) for c in range(K)])
    for a, b in E:
        for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
    for c in range(K): s.add_clause([-X(A, c), X(Bi, c)]); s.add_clause([X(A, c), -X(Bi, c)])
    if tri:
        for k, v in enumerate(tri): s.add_clause([X(v, k)])
    s.conf_budget(100_000_000)
    r = s.solve_limited(); st = s.accum_stats()
    if r is not True:
        save("FORCED_APART" if r is False else "hard")
        print(f"  iter {it}: n={n}: c(a) = c(b) is {'UNSAT -- a, b FORCED APART (verify with scripts/verify_gadget.py!)' if r is False else 'beyond budget'}   [{time.time()-t0:.0f}s]", flush=True)
        break
    m = s.get_model(); s.delete()
    col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
    rb = sorted(((0.25 * len(nb) - WN * near(p), k) for k, (p, nb) in pool.items()
                 if len(nb) >= 5 and len({col[j] for j in nb}) == K), reverse=True)
    if it % 5 == 1 or not rb:
        print(f"  iter {it}: n={n} m={len(E)} conflicts {st.get('conflicts', 0)}; {len(rb)} rainbow points   [{time.time()-t0:.0f}s]", flush=True)
    if it % 5 == 0: save("checkpoint")
    if not rb:
        save("pool_exhausted"); print("  pool exhausted for this colouring"); break
    for _, k in rb[:R]:
        p = pool.pop(k)[0]; vkey[k] = len(V); V.append(p)
    for i in range(n, len(V)): feed(i)
    G = build_graph(V)
