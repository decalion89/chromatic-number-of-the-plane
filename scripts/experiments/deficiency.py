"""delta(G): how many extra pairs must be forbidden to destroy five colours.

The 803-graph forces no single distance -- it stays 5-colourable with 1 and any
one d forbidden, 3500 extra edges included -- but it forces a PAIR of them:

    chi(five_247_c; {1, 0.542573, 0.209240}) > 5

so every 5-colouring of it has a monochromatic pair at one of those two
distances.  That is the project's first unconditional five-colour forcing
statement, and it suggests the metric that has been missing all along:

    delta(G) = fewest non-unit pairs whose forbidding kills 5-colourability

delta = 0 is exactly chi(G) >= 6.  Unlike mu, which has read 2 everywhere for
every graph here, delta is a number with room to move, it is monotone under
adding structure, and it can be driven down.

Computed the way certificates should be: give every extra pair a selector,
assume them all, and read the UNSAT core, which is a sufficient set in one
solve rather than thousands.  Re-assuming just the core shrinks it further, and
a greedy pass over what survives makes it minimal.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
POOL = int(sys.argv[2]) if len(sys.argv) > 2 else 40
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
print(f"{NAME} n={n} unit edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

buckets = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if dd > 16.0 or abs(dd - 1.0) < 1e-9: continue
        buckets[round(dd, 7)].append((i, j))
classes = sorted(buckets.items(), key=lambda kv: -len(kv[1]))[:POOL]
extra = []
for d2, prs in classes:
    for i, j in prs: extra.append((i, j, d2))
print(f"  pool: {len(classes)} classes, {len(extra)} candidate pairs   "
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
print(f"  all {len(extra)} forbidden: 5-colourable = {s.solve(assumptions=ass)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

core = set(s.get_core() or [])
sizes = [len(core)]
for _ in range(30):
    a2 = sorted(core)
    if s.solve(assumptions=a2): break
    c2 = set(s.get_core() or [])
    if len(c2) >= len(core): break
    core = c2; sizes.append(len(core))
    print(f"    core -> {len(core)}   [{time.time()-t0:.0f}s]", flush=True)
cur = sorted(core)
print(f"  core sizes {sizes}; greedy pass over {len(cur)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
i = 0
while i < len(cur):
    trial = cur[:i] + cur[i+1:]
    if trial and not s.solve(assumptions=trial):
        c2 = set(s.get_core() or trial)
        cur = [a for a in trial if a in c2] or trial
        i = 0
    else:
        i += 1
    if i and i % 25 == 0:
        print(f"    {len(cur)} left, probe {i}   [{time.time()-t0:.0f}s]", flush=True)
idx = [a - (n * K + 1) for a in cur]
by = defaultdict(int)
for t in idx: by[extra[t][2]] += 1
print(f"\n  delta <= {len(cur)} extra pairs   [{time.time()-t0:.0f}s]", flush=True)
print(f"  by distance: {dict(sorted(by.items(), key=lambda kv: -kv[1]))}", flush=True)
json.dump({"graph": NAME, "delta_upper_bound": len(cur),
           "pairs": [[extra[t][0], extra[t][1], extra[t][2]] for t in idx]},
          open(f"{ROOT}/data/deficiency_{NAME}", "w"))
