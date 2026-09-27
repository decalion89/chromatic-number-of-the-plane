"""Incremental colouring-guided growth: one persistent solver, learned clauses kept.

Adding vertices and edges only adds clauses, so a single incremental CaDiCaL
keeps everything it has learned from one iteration to the next (the older
scripts rebuilt graph and solver every time).  Each iteration: solve, read the
colouring, insert the R richest pool points whose neighbourhoods show all five
colours (the colouring cannot be extended to them), add only the new clauses.

MODE plain : grow until the graph itself refuses five colours.
MODE apart : clauses c(a) = c(b) kept permanently -- UNSAT means a, b forced apart.
MODE same  : clause c(a) != c(b) kept permanently -- UNSAT means a, b forced alike.
Every UNSAT is only a claim until verify6 / verify_gadget re-check it.

usage: grow_inc.py <in.json> <out.json> [R] [near_weight]   (env MODE, SOLVER, BUDGET)
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, os
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
MODE = os.environ.get("MODE", "plain"); SOLVER = os.environ.get("SOLVER", "cadical195")
BUDGET = int(os.environ.get("BUDGET", "200000000"))
d = json.load(open(IN))
F = Field(tuple(d["field_generators"]))
V = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
A, Bi = d.get("A"), d.get("B")
key = lambda p: (round(p.fx, 9), round(p.fy, 9))
vkey = {key(p): i for i, p in enumerate(V)}
G = build_graph(V)
E = list(G.edges())
U = {}
for a, b in E:
    for u in (V[b] - V[a], V[a] - V[b]): U.setdefault(key(u), u)
U = list(U.values())
adj = [set() for _ in range(len(V))]
for a, b in E: adj[a].add(b); adj[b].add(a)
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
    json.dump({"field_generators": list(F.gens), "A": A, "B": Bi, "status": tag, "mode": MODE,
               "points": [[[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]] for q in V]},
              open(OUT, "w"))
X = lambda v, c: 1 + v * K + c
s = Solver(name=SOLVER)
def add_vertex(v):
    s.add_clause([X(v, c) for c in range(K)])
def add_edge(a, b):
    for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
for v in range(len(V)): add_vertex(v)
for a, b in E: add_edge(a, b)
tri = next(((a, b, c) for a, b in E for c in sorted(adj[a] & adj[b])), None)
if tri:
    for k, v in enumerate(tri): s.add_clause([X(v, k)])
if MODE == "apart":
    for c in range(K): s.add_clause([-X(A, c), X(Bi, c)]); s.add_clause([X(A, c), -X(Bi, c)])
elif MODE == "same":
    for c in range(K): s.add_clause([-X(A, c), -X(Bi, c)])
if A is not None:
    ax, ay, bx, by = V[A].fx, V[A].fy, V[Bi].fx, V[Bi].fy
    near = lambda p: min(((p.fx - ax) ** 2 + (p.fy - ay) ** 2) ** 0.5, ((p.fx - bx) ** 2 + (p.fy - by) ** 2) ** 0.5)
else:
    cx = sum(p.fx for p in V) / len(V); cy = sum(p.fy for p in V) / len(V)
    near = lambda p: ((p.fx - cx) ** 2 + (p.fy - cy) ** 2) ** 0.5
m_edges = len(E)
print(f"{IN}: n={len(V)} m={m_edges}, {len(U)} directions, pool {len(pool)}; mode {MODE}, solver {SOLVER}"
      + (f", pair {A},{Bi} at d^2 = {V[A].dist2(V[Bi])}" if A is not None else "") + f"   [{time.time()-t0:.0f}s]", flush=True)
last = 0
for it in range(1, 1000000):
    n = len(V)
    s.conf_budget(BUDGET)
    r = s.solve_limited(); st = s.accum_stats(); conf = st.get("conflicts", 0) - last; last = st.get("conflicts", 0)
    if r is not True:
        tag = ("UNSAT_" + MODE) if r is False else "hard"
        save(tag)
        print(f"  iter {it}: n={n} m={m_edges}: {'UNSAT (' + MODE + ') -- a claim until verified' if r is False else 'beyond budget'}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        break
    mdl = s.get_model()
    col = [next(c for c in range(K) if mdl[X(v, c) - 1] > 0) for v in range(n)]
    rb = sorted(((0.25 * len(nb) - WN * near(p), k) for k, (p, nb) in pool.items()
                 if len(nb) >= 5 and len({col[j] for j in nb}) == K), reverse=True)
    if it % 5 == 1 or not rb:
        print(f"  iter {it}: n={n} m={m_edges} conflicts {conf}; {len(rb)} rainbow points   [{time.time()-t0:.0f}s]", flush=True)
    if it % 5 == 0: save("checkpoint")
    if not rb:
        save("pool_exhausted"); print("  pool exhausted for this colouring", flush=True); break
    for _, k in rb[:R]:
        p, nb = pool.pop(k); v = len(V); vkey[k] = v; V.append(p); adj.append(set())
        add_vertex(v)
        # neighbours: every existing vertex at a unit step along the module's directions
        for u in U:
            j = vkey.get(key(p + u))
            if j is not None and j != v and j not in adj[v]:
                adj[v].add(j); adj[j].add(v); add_edge(v, j); m_edges += 1
    for i in range(n, len(V)): feed(i)
