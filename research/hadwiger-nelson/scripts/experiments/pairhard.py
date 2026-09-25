"""Is the pair itself constrained, or is the whole graph just hard?

Solve the saved graph three ways with the same pinned triangle and compare the
solver's work: plain 5-colouring, with c(a) = c(b), with c(a) != c(b).  A pair
close to forced apart shows up as the equal case alone growing expensive.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
K = 5
d = json.load(open(sys.argv[1]))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
a, b = d["A"], d["B"]
g = build_graph(P); n = g.n; E = list(g.edges())
adj = [set() for _ in range(n)]
for x, y in E: adj[x].add(y); adj[y].add(x)
tri = next((x, y, z) for x, y in E for z in adj[x] & adj[y])
X = lambda v, c: 1 + v * K + c
print(f"{sys.argv[1]}: n={n} m={len(E)}; deg a={len(adj[a])} deg b={len(adj[b])}", flush=True)
for label, extra in (("plain", []),
                     ("c(a) = c(b)", [cl for c in range(K) for cl in ([-X(a, c), X(b, c)], [X(a, c), -X(b, c)])]),
                     ("c(a) != c(b)", [[-X(a, c), -X(b, c)] for c in range(K)])):
    for seed in range(3):
        s = Solver(name="cd19")
        for v in range(n): s.add_clause([X(v, c) for c in range(K)])
        for x, y in E:
            for c in range(K): s.add_clause([-X(x, c), -X(y, c)])
        for cl in extra: s.add_clause(cl)
        for k, v in enumerate(tri): s.add_clause([X(v, k)])
        # perturb: pin one more vertex to a colour different from its pinned neighbours' by seed
        t0 = time.time(); s.conf_budget(20_000_000); r = s.solve_limited()
        print(f"  {label:<13s} seed {seed}: {r}  conflicts {s.accum_stats().get('conflicts')}  [{time.time()-t0:.1f}s]", flush=True)
        s.delete()
        break
