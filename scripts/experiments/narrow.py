"""The narrowest pigeonhole-free disjunction: few PAIRS, but only two distances.

Asking for the fewest individual pairs outright gives ten, and they are
pigeonhole: six vertices, five unit edges between them, forbid the rest and you
have K6.  Restricting every candidate pair to ONE distance class makes that
impossible -- a two-distance set in the plane has at most five points -- but a
whole class forbidden still leaves the 803-graph 5-colourable, so no subset can
work either.

Two classes is the sweet spot.  Forbidding both of

    d^2 = 3/2 - sqrt33/6   and   d^2 = 7/6 - sqrt33/6

kills five colours, and the resulting three-distance graph was checked for a
6-clique and has none -- so neither does any subgraph of it.  Every subset of
those 4728 pairs is therefore pigeonhole-free by inheritance, and the minimal
one is a genuine narrow disjunction:

    in every 5-colouring of the 803-graph, one of THESE pairs is monochromatic.

Selectors and an UNSAT core find it in one solve rather than thousands, and a
greedy pass makes it minimal.  Then the only question left is whether the
survivors share any structure a rotation could consume -- a common vertex, or
an orbit.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, itertools
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = "five_247_c.json"
D2 = [0.5425729, 0.2092396]
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
extra = []
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        for t in D2:
            if abs(dd - t) < 1e-7:
                extra.append((i, j, t)); break
print(f"{NAME} n={n} unit={len(E)} candidate pairs={len(extra)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

X = lambda v, c: 1 + v * K + c
S = lambda t: n * K + 1 + t
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])
for t, (i, j, _) in enumerate(extra):
    for c in range(K):
        cnf.append([-S(t), -X(i, c), -X(j, c)])
s = Solver(name="cd19", bootstrap_with=cnf)
ass = [S(t) for t in range(len(extra))]
print(f"  all forbidden: 5-colourable = {s.solve(assumptions=ass)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
core = set(s.get_core() or [])
print(f"  first core: {len(core)}   [{time.time()-t0:.0f}s]", flush=True)
for _ in range(40):
    a2 = sorted(core)
    if s.solve(assumptions=a2): break
    c2 = set(s.get_core() or [])
    if len(c2) >= len(core): break
    core = c2
    print(f"    core -> {len(core)}   [{time.time()-t0:.0f}s]", flush=True)
cur = sorted(core); i = 0
while i < len(cur):
    trial = cur[:i] + cur[i+1:]
    if trial and not s.solve(assumptions=trial):
        c2 = set(s.get_core() or trial)
        cur = [a for a in trial if a in c2] or trial
        i = 0
    else:
        i += 1
idx = [a - (n * K + 1) for a in cur]
pairs = [extra[t] for t in idx]
print(f"\n  MINIMAL DISJUNCTION: {len(pairs)} pairs   [{time.time()-t0:.0f}s]",
      flush=True)
vs = Counter()
for i0, j0, t in pairs: vs[i0] += 1; vs[j0] += 1
print(f"  distances: {Counter(round(t,7) for _,_,t in pairs)}", flush=True)
print(f"  {len(vs)} distinct vertices; busiest: {vs.most_common(6)}", flush=True)
for i0, j0, t in pairs:
    print(f"     ({i0},{j0})  d^2={t:.7f}", flush=True)
# does any single vertex meet them all?  that would be a hub disjunction
full = [v for v, c in vs.items() if c == len(pairs)]
print(f"  vertices meeting every pair (a hub): {full}", flush=True)
json.dump({"graph": NAME, "pairs": [[a, b, t] for a, b, t in pairs]},
          open(f"{ROOT}/data/narrow_disjunction.json", "w"))
