"""Give the 5-chromatic carrier the symmetry its glues need.

Overlap is the resource -- with the copies nearly disjoint, colouring the
second by sigma . c satisfies every cross edge and nothing is forced.  And the
measured deficit is exactly there: Sa's best glue reuses 56% of it and, after
a level of thickening, 85%, while the 951-vertex 5-chromatic graph tops out at
39%.  The reason is plain.  Sa is a D6 orbit; the spindled graphs are not,
because the spindle rotation is about a vertex and breaks the symmetry.

So restore it: take the orbit of the 5-chromatic graph under the 12-element
group Sa already carries.  Sa itself is fixed by that group, so the orbit costs
far less than twelve times the size, and whatever symmetry it gains is exactly
what the high-overlap glues need.
"""
import sys, json, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

SRC = sys.argv[1]
K = 5
d = json.load(open(SRC))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
print(f"{SRC}: {len(P)} points", flush=True)
rot60 = _rot60(F)
seen, S = set(), []
for p in P:
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            if q not in seen: seen.add(q); S.append(q)
            q = rot60(q)
t0 = time.time()
g = build_graph(S); n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"D6 orbit: n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]", flush=True)

X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf)
if not s.solve():
    print("  *** THE ORBIT REFUSES FIVE COLOURS ***", flush=True); sys.exit()
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
r = random.Random(4)
for _ in range(60):
    if len(cols) >= 20: break
    s.set_phases([(1 if r.random() < 0.2 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): continue
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
print(f"  free@5={100.0*free/n:.2f}%  forced-pair candidates {len(cand)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
conf = 0
for a, b in cand[:4000]:
    s2 = Solver(name="m22", bootstrap_with=cnf)
    diff = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
    if not diff:
        conf += 1
        dd = (S[a] - S[b]).norm2()
        print(f"    FORCED v{a} v{b}  d^2={dd}  <<<<<<<<<<<<", flush=True)
print(f"  confirmed forced pairs at five: {conf}", flush=True)

SS = set(S); best = []
for name, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3)), ("rot180", Fr(1, 4))):
    r0 = rotation_joining(d2, F)
    for c in S[:400]:
        rot = r0.about(c)
        ov = sum(1 for p in S if rot(p) in SS)
        if 1 < ov < n: best.append(ov)
best.sort(reverse=True)
print(f"  best two-copy glue overlaps: {best[:6]} of {n} "
      f"({100.0*best[0]/n:.0f}%)   [{time.time()-t0:.0f}s]", flush=True)
