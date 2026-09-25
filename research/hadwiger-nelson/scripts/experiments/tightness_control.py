"""The control the tightness numbers need: the same six copies, NOT touching.

Six copies of the 803 coupled on a hexagon need a median of 207 370 conflicts to
5-colour, against 33 for one copy.  But the coupled graph is five times larger,
and conflicts grow with size for reasons that have nothing to do with how
constrained a graph is.  The right baseline is the same six copies placed far
apart -- same vertices, same edges within copies, no coupling between them -- so
whatever difference remains is the coupling and nothing else.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random, statistics
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
pts = []
for i in range(6):
    off = F.rational(100 * i)
    pts += [Point(p.x + off, p.y) for p in P]
G = build_graph(pts); n = G.n; E = list(G.edges())
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for a, b in E:
    for c in range(K): cnf.append([-X(a, c), -X(b, c)])
rng = random.Random(0); out = []
for _ in range(3):
    s = Solver(name="cd19", bootstrap_with=cnf)
    v = rng.randrange(n); c = rng.randrange(K)
    ok = s.solve(assumptions=[X(v, c)]); st = s.accum_stats(); s.delete()
    out.append(st.get("conflicts", 0))
print(f"  six disjoint copies: n={n} edges={len(E)} conflicts {out} "
      f"median {statistics.median(out)}   [{time.time()-t0:.0f}s]", flush=True)
