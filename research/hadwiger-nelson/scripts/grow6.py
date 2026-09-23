"""Colouring-guided growth: insert the points each colouring cannot extend to.

Counterexample-guided search for a non-5-colourable unit-distance graph.  Hold a
graph G and a 5-colouring c of it.  A candidate point p (a unit step from some
vertex, along a direction G already uses) whose G-neighbours show all five
colours under c cannot be coloured: adding p kills c and every other colouring
in which that neighbourhood is a rainbow.  Add the richest such points, re-solve,
repeat.  UNSAT would be a 6-chromatic candidate (verify before believing);
running dry means the current colouring extends point by point to the pool.
"""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
R = int(sys.argv[2]) if len(sys.argv) > 2 else 40
ITERS = int(sys.argv[3]) if len(sys.argv) > 3 else 400
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
G = build_graph(V)
U = {}
for a, b in G.edges():
    for u in (V[b] - V[a], V[a] - V[b]): U.setdefault(key(u), u)
U = list(U.values())
print(f"{NAME}: n={G.n} m={G.m}, {len(U)} edge directions   [{time.time()-t0:.0f}s]", flush=True)
vkey = {key(p): i for i, p in enumerate(V)}
pool = {}
def feed(i):
    for u in U:
        p = V[i] + u; k = key(p)
        if k in vkey: continue
        e = pool.get(k)
        if e is None: pool[k] = [p, {i}]
        else: e[1].add(i)
for i in range(len(V)): feed(i)
print(f"  pool {len(pool)} points, {sum(1 for e in pool.values() if len(e[1]) >= 5)} with 5+ neighbours"
      f"   [{time.time()-t0:.0f}s]", flush=True)
for it in range(1, ITERS + 1):
    n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    s = Solver(name="cd19")
    for v in range(n): s.add_clause([X(v, c) for c in range(K)])
    for a, b in E:
        for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
    s.conf_budget(30_000_000)
    r = s.solve_limited(); st = s.accum_stats()
    if r is None:
        print(f"  iter {it}: n={n} m={len(E)}: solver budget exhausted -- hard instance, saving", flush=True)
        r = "unknown"
    if r is False or r == "unknown":
        json.dump({"field_generators": list(F.gens), "status": str(r),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]] for q in V]},
                  open(f"{ROOT}/data/grow6_{'unsat' if r is False else 'hard'}_{n}.json", "w"))
        if r is False:
            print(f"  *** iter {it}: n={n} NOT 5-colourable -- VERIFY INDEPENDENTLY ***", flush=True)
        break
    m = s.get_model(); s.delete()
    col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
    rb = [(len(nb), k) for k, (p, nb) in pool.items()
          if len(nb) >= 5 and len({col[j] for j in nb}) == K]
    if it % 5 == 1 or not rb:
        print(f"  iter {it}: n={n} m={len(E)} conflicts {st.get('conflicts', 0)}; "
              f"{len(rb)} rainbow candidates, richest {sorted(rb, reverse=True)[:1]}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
    if not rb:
        print("  the colouring extends to every pool point one at a time -- pool exhausted")
        break
    rb.sort(reverse=True)
    for _, k in rb[:R]:
        p = pool.pop(k)[0]; vkey[k] = len(V); V.append(p)
    for i in range(n, len(V)): feed(i)
    G = build_graph(V)
