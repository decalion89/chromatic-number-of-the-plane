"""The D6 orbit of the 951-vertex graph: slack, forcing, and glue overlap.

Two corrections over the first attempt.  Cadical settles the orbit's
5-colourability in four seconds where minisat ground for a quarter of an hour,
so the solving is done with cadical; and the diversity trick that makes the
forced-pair filter work -- randomising every decision polarity -- is what was
hurting it, because on a large loose instance it forces the solver to fight its
own heuristic.  Randomising a tenth of the polarities gives different
colourings without that cost.

The point of symmetrising is overlap.  Sa's best glue reuses 56% of it because
Sa is a D6 orbit; the spindled 5-chromatic graphs top out at 39% because the
spindle rotation is about a vertex and breaks the symmetry.  Restoring the
symmetry should restore the overlap, and overlap is what defeats the
sigma-argument.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

SRC = sys.argv[1]; K = int(sys.argv[2]); NCOL = int(sys.argv[3])
d = json.load(open(SRC))
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
g = build_graph(S); n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"{SRC}: {len(P)} -> D6 orbit n={n} m={m} deg={2.0*m/n:.2f}", flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
t0 = time.time()
s = Solver(name="cd19", bootstrap_with=cnf)
if not s.solve():
    print(f"  *** THE ORBIT REFUSES {K} COLOURS ***   [{time.time()-t0:.0f}s]", flush=True)
    sys.exit()
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
r = random.Random(4)
while len(cols) < NCOL:
    # set_phases steers minisat -- that is what makes the filter work on Sa
    # and on G -- but cadical ignores it, and twenty "different" colourings
    # came back identical, leaving 4.4 million pairs, which is no filter at
    # all.  Blocking a random sample of the last solution forces a genuinely
    # different one whatever the solver does with hints.
    last = cols[-1]
    samp = r.sample(range(n), 40)
    s.add_clause([-X(v, last[v]) for v in samp])
    if not s.solve(): break
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
free = min(sum(1 for v in range(n) if len({col[u] for u in g.adj[v]} | {col[v]}) < K)
           for col in cols)
buck = defaultdict(list)
for v in range(n):
    buck[tuple(c[v] for c in cols)].append(v)
cand = [(a, b) for vs in buck.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i+1:]]
print(f"  {len(cols)} colourings, free@{K}={100.0*free/n:.2f}%, "
      f"{len(cand)} forced-pair candidates   [{time.time()-t0:.0f}s]", flush=True)
conf = 0
for a, b in cand[:3000]:
    s2 = Solver(name="cd19", bootstrap_with=cnf)
    diff = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
    if not diff:
        conf += 1
        print(f"    FORCED v{a} v{b}  d^2={(S[a]-S[b]).norm2()}  <<<<<<<<", flush=True)
print(f"  confirmed forced pairs at {K}: {conf}   [{time.time()-t0:.0f}s]", flush=True)
SS = set(S); best = []
for name, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3)), ("rot180", Fr(1, 4))):
    r0 = rotation_joining(d2, F)
    for c in S[:300]:
        rot = r0.about(c)
        ov = sum(1 for p in S if rot(p) in SS)
        if 1 < ov < n: best.append(ov)
best.sort(reverse=True)
print(f"  best glue overlaps: {best[:6]} of {n} ({100.0*best[0]/n:.0f}%)"
      f"   [{time.time()-t0:.0f}s]", flush=True)
