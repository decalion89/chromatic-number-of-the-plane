"""Test the only pairs a 5-colour forcing can ever hold: u - v in 5M.

The admissible coset colourings of the tight graph span the whole dual of its
edge module, so their common kernel is exactly 5M: every pair whose difference
is NOT five times a module element is split by some coset colouring, and no
growth along these directions can ever force it together.  What is left is
distance 5 along a unit direction -- Exoo-Ismailescu's distance, closed by their
rotation (49 + 3 sqrt-11)/50, which carries 5 in its denominator and so escapes
every coset colouring of the union.  One solver call per surviving pair.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
ROOT = HN_DIR
NAME = sys.argv[1]; t0 = time.time(); K = 5
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n; E = list(G.edges())
cand = [(i, j) for i, j, T in json.load(open(f"{ROOT}/data/coset_unsplit_{NAME}"))["unsplittable_by_cosets"]]
print(f"{NAME}: n={n}; {len(cand)} coset-unsplittable pairs to settle", flush=True)
X = lambda v, c: 1 + v * K + c
s = Solver(name="cd19")
for v in range(n): s.add_clause([X(v, c) for c in range(K)])
for a, b in E:
    for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
nsel = n * K + 1; forced = []; calls = 0
while cand:
    i, j = cand[0]; sel = nsel; nsel += 1
    for c in range(K): s.add_clause([-sel, -X(i, c), -X(j, c)])
    r = s.solve(assumptions=[sel]); calls += 1
    st = s.accum_stats()
    m = s.get_model() if r else None
    s.add_clause([-sel])
    if not r:
        forced.append((i, j))
        print(f"  *** ({i},{j}) FORCED TOGETHER AT FIVE, distance 5 -- VERIFY ***   [{time.time()-t0:.0f}s]", flush=True)
        json.dump({"graph": NAME, "forced_same_5M": forced}, open(f"{ROOT}/data/forced5M_{NAME}", "w"))
        cand = cand[1:]; continue
    col = [next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)]
    before = len(cand)
    cand = [(a, b) for a, b in cand if col[a] == col[b]]
    print(f"    call {calls}: split ({i},{j}) and {before - len(cand) - 1} more; {len(cand)} left  "
          f"(conflicts so far {st.get('conflicts', 0)})   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  forced-together 5M pairs: {len(forced)}   [{time.time()-t0:.0f}s]", flush=True)
