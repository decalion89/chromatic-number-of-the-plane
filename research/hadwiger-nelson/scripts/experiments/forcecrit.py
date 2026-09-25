"""Is EVERY vertex needed to force the pair, or only the last one in the order?

The prefix search returned 802 of 802 for every genuinely forced pair at the
degree-4 ceiling point, but a prefix is one nested family among many: it shows
that the BFS-last vertex is needed, not that all of them are.  The theorem
wanted is stronger and is one pass of single deletions:

    for every w, is (u, v) still forced apart in H - w ?

If the answer is no for every w, then H is vertex-critical FOR THIS FORCING,
and no subset of it certifies the pair -- the minimal certificate is the whole
near-critical graph, exactly.  Each deletion that breaks the forcing makes the
instance satisfiable, which is the cheap direction, so the pass is affordable.
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

ROOT = HN_DIR
t0 = time.time(); K = 4
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
print(f"five_247_c n={n} m={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

U, V, PIVOT = 461, 561, 315          # d^2 = 1/3, the disjunctive-spindle distance
X = lambda w, c: 1 + w * K + c
def forced_without(drop):
    ks = set(range(n)) - drop
    cnf = [[X(w, c) for c in range(K)] for w in ks]
    for x, y in E:
        if x in ks and y in ks:
            for c in range(K):
                cnf.append([-X(x, c), -X(y, c)])
    cnf.append([X(U, 0)]); cnf.append([X(V, 0)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    r = s.solve(); s.delete()
    return not r

print(f"pair ({U},{V}) at d^2 = 1/3, pivot {PIVOT}", flush=True)
print(f"  forced in H: {forced_without({PIVOT})}   [{time.time()-t0:.0f}s]", flush=True)
needed = []; spare = []
for w in range(n):
    if w in (U, V, PIVOT): continue
    if forced_without({PIVOT, w}):
        spare.append(w)
    else:
        needed.append(w)
    if (len(needed) + len(spare)) % 100 == 0:
        print(f"    {len(needed)+len(spare)} tested: {len(needed)} needed, "
              f"{len(spare)} spare   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  vertices whose deletion BREAKS the forcing: {len(needed)} of {n-3}",
      flush=True)
print(f"  vertices that are spare: {len(spare)}   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"pair": [U, V], "pivot": PIVOT, "needed": len(needed),
           "spare": spare}, open(f"{ROOT}/data/forcing_criticality.json", "w"))
