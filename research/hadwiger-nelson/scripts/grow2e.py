"""Targeted growth for the Exoo-Ismailescu gadget: force u, u + 2e APART.

The 803-graph's module passes the apart gate in exactly one direction among
those tested: e = (-3/16, sqrt247/16), the single sqrt247 edge, which lies in
8M.  Take that edge a -- a + e, adjoin b = a + 2e, and repeat: ask the solver
for a 5-colouring with c(a) = c(b); add the pool points whose neighbourhoods
it shows all five colours on; rebuild.  UNSAT would be a unit-distance graph
splitting a distance-2 pair at five colours -- with Exoo-Ismailescu, a proof
of chi >= 6.  Verify before believing anything.
"""
import sys, time, json, math
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
R = int(sys.argv[1]) if len(sys.argv) > 1 else 20
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(V)
e26 = Point(F.rational(Fr(-3, 16)), F.rational(Fr(1, 16)) * F.sqrt(247))
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
vkey = {key(p): i for i, p in enumerate(V)}
A = next(i for i, p in enumerate(V) if key(p + e26) in vkey or key(p - e26) in vkey)
if key(V[A] + e26) not in vkey: e26 = Point(-e26.x, -e26.y)
Bp = V[A] + e26 + e26
if key(Bp) not in vkey:
    vkey[key(Bp)] = len(V); V.append(Bp)
Bi = vkey[key(Bp)]
assert V[A].dist2(V[Bi]) == F.rational(4)
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
print(f"803 + b: n={G.n} m={G.m}; pair {A},{Bi} at distance 2 along the sqrt247 edge; pool {len(pool)}   [{time.time()-t0:.0f}s]", flush=True)
for it in range(1, 100000):
    n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    s = Solver(name="cd19")
    for v in range(n): s.add_clause([X(v, c) for c in range(K)])
    for a, b in E:
        for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
    for c in range(K): s.add_clause([-X(A, c), X(Bi, c)]); s.add_clause([X(A, c), -X(Bi, c)])
    s.conf_budget(50_000_000)
    r = s.solve_limited(); st = s.accum_stats()
    if r is not True:
        tag = "FORCED_APART" if r is False else "hard"
        json.dump({"field_generators": list(F.gens), "A": A, "B": Bi, "status": tag,
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in V]},
                  open(f"{ROOT}/data/grow2e_{tag}.json", "w"))
        print(f"  iter {it}: n={n}: c(a) = c(b) is {'UNSAT -- a, b FORCED APART (verify!)' if r is False else 'beyond budget'}", flush=True)
        break
    m = s.get_model(); s.delete()
    col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
    rb = sorted(((len(nb), k) for k, (p, nb) in pool.items() if len(nb) >= 5 and len({col[j] for j in nb}) == K), reverse=True)
    if it % 5 == 1 or not rb:
        print(f"  iter {it}: n={n} m={len(E)} conflicts {st.get('conflicts', 0)}; {len(rb)} rainbow points   [{time.time()-t0:.0f}s]", flush=True)
    if not rb:
        print("  pool exhausted for this colouring"); break
    for _, k in rb[:R]:
        p = pool.pop(k)[0]; vkey[k] = len(V); V.append(p)
    for i in range(n, len(V)): feed(i)
    G = build_graph(V)
