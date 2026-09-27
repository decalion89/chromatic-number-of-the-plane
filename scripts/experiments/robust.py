"""How robust is E-I's forcing?  Minimum number of H's two-edges that must be ALIKE in a proper
unit-edge 5-colouring of H with c(A) != c(B) (MaxSAT, RC2).  If it is m, a unit-distance graph
only has to forbid m simultaneous alike two-edges -- not force every one apart."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, sys, time
sys.path = [p for p in sys.path if "scratchpad" not in p]   # an old six.py there shadows the six package
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.formula import WCNF
from pysat.examples.rc2 import RC2
d = json.load(open(HN_DIR + "/data/ei_H214.json")); F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
A, B = d["A"], d["B"]; n = len(P)
g = build_graph(P); E = list(g.edges())
T = [(i, j) for i in range(n) for j in range(i + 1, n) if P[i].dist2(P[j]) == 4]
X = lambda v, c: 1 + 5 * v + c
w = WCNF(); top = 5 * n
for v in range(n): w.append([X(v, c) for c in range(5)])
for a, b in E:
    for c in range(5): w.append([-X(a, c), -X(b, c)])
for c in range(5): w.append([-X(A, c), -X(B, c)])
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(b); adj[b].add(a)
tri = next(((a, b, c) for a, b in E for c in sorted(adj[a] & adj[b])), None)
for k, v in enumerate(tri): w.append([X(v, k)])
sel = []
for i, j in T:
    top += 1; s = top; sel.append(s)
    for c in range(5): w.append([-s, -X(i, c), -X(j, c)])     # s -> two-edge (i, j) not alike
    w.append([s], weight=1)
t0 = time.time()
with RC2(w, solver="g4") as rc:
    m = rc.compute()
    print(f"H: {n} points, {len(E)} unit edges, {len(T)} two-edges; min #alike two-edges with c(A) != c(B): {rc.cost}   [{time.time()-t0:.1f}s]")
    alike = [T[k] for k, s in enumerate(sel) if m[s - 1] < 0]
    print("  alike two-edges in an optimum:", alike[:12])
