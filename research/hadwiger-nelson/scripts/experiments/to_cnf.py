"""Write the 5-colouring CNF of a saved growth instance (edges along the stored unit vectors,
MODE constraint included, one triangle pinned) as DIMACS, for an independent solver + DRAT."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
d = json.load(open(sys.argv[1])); out = sys.argv[2]; K = 5
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]
key = lambda x, y: (round(x, 9), round(y, 9))
vk = {key(p.fx, p.fy): i for i, p in enumerate(V)}
Ufx = [(u.fx, u.fy) for u in U]
E = set()
for i, p in enumerate(V):
    for ux, uy in Ufx:
        j = vk.get(key(p.fx + ux, p.fy + uy))
        if j is not None and j != i: E.add((min(i, j), max(i, j)))
E = sorted(E); n = len(V)
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(b); adj[b].add(a)
tri = next(((a, b, c) for a, b in E for c in sorted(adj[a] & adj[b])), None)
X = lambda v, c: 1 + v * K + c
cl = [[X(v, c) for c in range(K)] for v in range(n)]
for a, b in E:
    for c in range(K): cl.append([-X(a, c), -X(b, c)])
mode = d.get("mode", "plain"); A, B = d.get("A"), d.get("B")
if mode == "apart":
    for c in range(K): cl.append([-X(A, c), X(B, c)]); cl.append([X(A, c), -X(B, c)])
elif mode == "same":
    for c in range(K): cl.append([-X(A, c), -X(B, c)])
for k, v in enumerate(tri): cl.append([X(v, k)])
with open(out, "w") as f:
    f.write(f"c {sys.argv[1]} n={n} m={len(E)} mode={mode} triangle={tri}\n")
    f.write(f"p cnf {n * K} {len(cl)}\n")
    for c in cl: f.write(" ".join(map(str, c)) + " 0\n")
print(f"{out}: n={n} m={len(E)} mode={mode}; {n*K} vars, {len(cl)} clauses; triangle {tri}")
