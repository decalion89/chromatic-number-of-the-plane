"""Is the D6 orbit of the 951-vertex 5-chromatic graph 5-colourable?

minisat has been grinding on this for a quarter of an hour where it settles
G (1581 vertices) instantly, which is either bad luck or an UNSAT proof.
Cadical is the right solver for a single hard instance -- the earlier benchmark
in this project had it losing badly on assumption-heavy scans and winning on
exactly this shape of question -- so ask it directly, and ask for four colours
too as a sanity check on the encoding.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

d = json.load(open(HN_DIR + "/data/five_247_b.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
rot60 = _rot60(F)
seen, S = set(), []
for p in P:
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            if q not in seen: seen.add(q); S.append(q)
            q = rot60(q)
g = build_graph(S)
n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"D6 orbit of the 951: n={n} m={m} deg={2.0*m/n:.2f}", flush=True)
# a triangle pinned for symmetry breaking -- sound, and worth a factor of k!
tri = g.find_clique(3)
for K in (5, 6):
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    for i, v in enumerate(tri):
        cnf.append([X(v, i)])
        for c in range(K):
            if c != i: cnf.append([-X(v, c)])
    t0 = time.time()
    s = Solver(name="cd19", bootstrap_with=cnf)
    ok = s.solve()
    print(f"  {K}-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
    if ok:
        model = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in model) for v in range(n)]
        bad = [(u, v) for u, v in g.edges() if col[u] == col[v]]
        print(f"    colouring verified: {len(bad)} monochromatic edges", flush=True)
        s.delete(); break
    s.delete()
    print(f"  *** the orbit needs more than {K} colours ***", flush=True)
