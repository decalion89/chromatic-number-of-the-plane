"""Targeted growth: kill every 5-colouring that splits one pair in 5M.

In an integral multiquadratic graph the coset colourings never split a pair
whose difference lies in 5M, and growth along the graph's own directions never
kills a coset colouring.  So hold a pair A, B = A + 5e and repeat: ask the
solver for a 5-colouring with c(A) != c(B); add the pool points (one unit step
along an existing direction) whose neighbourhoods that colouring shows all five
colours on; rebuild.  If the question ever goes UNSAT, A and B are forced
together, and Exoo-Ismailescu's rotation (49 + 3 sqrt-11)/50 about A closes
the distance 5 into a candidate 6-chromatic graph -- verify before believing.
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
NAME = sys.argv[1]; R = int(sys.argv[2]) if len(sys.argv) > 2 else 20
ITERS = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
A, B = d["A"], d["B"]
assert V[A].dist2(V[B]) == F.rational(25)
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
G = build_graph(V)
U = {}
for a, b in G.edges():
    for u in (V[b] - V[a], V[a] - V[b]): U.setdefault(key(u), u)
U = list(U.values())
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
print(f"{NAME}: n={G.n} m={G.m}, {len(U)} directions, pool {len(pool)}; target pair {A},{B} at distance 5"
      f"   [{time.time()-t0:.0f}s]", flush=True)
hist = []
for it in range(1, ITERS + 1):
    n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    s = Solver(name="cd19")
    for v in range(n): s.add_clause([X(v, c) for c in range(K)])
    for a, b in E:
        for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
    for c in range(K): s.add_clause([-X(A, c), -X(B, c)])
    s.conf_budget(50_000_000)
    r = s.solve_limited(); st = s.accum_stats()
    if r is not True:
        tag = "FORCED" if r is False else "hard"
        json.dump({"field_generators": list(F.gens), "A": A, "B": B, "status": tag,
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]] for q in V]},
                  open(f"{ROOT}/data/grow5m_{tag}_{NAME}", "w"))
        print(f"  iter {it}: n={n} m={len(E)}: splitting A,B is {'UNSAT -- A, B FORCED TOGETHER (verify!)' if r is False else 'beyond the budget'}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        break
    m = s.get_model(); s.delete()
    col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
    rb = sorted(((len(nb), k) for k, (p, nb) in pool.items()
                 if len(nb) >= 5 and len({col[j] for j in nb}) == K), reverse=True)
    hist.append(st.get("conflicts", 0))
    if it % 5 == 1 or not rb:
        print(f"  iter {it}: n={n} m={len(E)} conflicts {hist[-1]}; {len(rb)} rainbow points"
              f"   [{time.time()-t0:.0f}s]", flush=True)
    if not rb:
        print("  no pool point is a rainbow in this splitting colouring: the split extends to the whole pool")
        break
    for _, k in rb[:R]:
        p = pool.pop(k)[0]; vkey[k] = len(V); V.append(p)
    for i in range(n, len(V)): feed(i)
    G = build_graph(V)
