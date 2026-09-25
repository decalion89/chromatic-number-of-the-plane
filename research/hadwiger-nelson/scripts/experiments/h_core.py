"""Which of H's 446 two-edges does Exoo-Ismailescu's forcing actually need?
Selector literal per two-edge, solve 'c(A) != c(B)' under all selectors, shrink the UNSAT core
by re-solving under the core, then try deleting each remaining two-edge (deletion-based MUS).
A small core = fewer places where a unit-distance gadget has to do the two-edge's job."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, sys, time
sys.path = [p for p in sys.path if "scratchpad" not in p]
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
d = json.load(open(HN_DIR + "/data/ei_H214.json")); F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
A, B = d["A"], d["B"]; n = len(P)
g = build_graph(P); E = list(g.edges())
T = [(i, j) for i in range(n) for j in range(i + 1, n) if P[i].dist2(P[j]) == 4]
X = lambda v, c: 1 + 5 * v + c
cl = [[X(v, c) for c in range(5)] for v in range(n)]
for a, b in E:
    for c in range(5): cl.append([-X(a, c), -X(b, c)])
for c in range(5): cl.append([-X(A, c), -X(B, c)])
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(b); adj[b].add(a)
tri = next(((a, b, c) for a, b in E for c in sorted(adj[a] & adj[b])), None)
for k, v in enumerate(tri): cl.append([X(v, k)])
sel = []
top = 5 * n
for i, j in T:
    top += 1; sel.append(top)
    for c in range(5): cl.append([-top, -X(i, c), -X(j, c)])
s = Solver(name="cadical195", bootstrap_with=cl)
t0 = time.time()
core = list(sel)
for rnd in range(6):
    r = s.solve(assumptions=core)
    assert r is False, "not UNSAT under the current selectors"
    new = s.get_core()
    print(f"  round {rnd}: core {len(new)} of {len(core)} two-edges   [{time.time()-t0:.0f}s]", flush=True)
    if len(new) >= len(core): break
    core = sorted(new)
# deletion-based minimisation
i = 0; kept = list(core)
while i < len(kept):
    trial = kept[:i] + kept[i + 1:]
    if s.solve(assumptions=trial) is False:
        c2 = s.get_core(); kept = [x for x in trial if x in set(c2)]
    else:
        i += 1
    print(f"  deletion: {len(kept)} two-edges left, position {i}   [{time.time()-t0:.0f}s]", flush=True)
idx = {v: k for k, v in enumerate(sel)}
mus = [T[idx[x]] for x in kept]
json.dump({"mus_two_edges": mus, "A": A, "B": B}, open(HN_DIR + "/data/ei_H_two_edge_mus.json", "w"))
print(f"MUS: {len(mus)} of {len(T)} two-edges suffice for E-I's forcing   [{time.time()-t0:.0f}s]")
