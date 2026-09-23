"""The densest object the project can build, carried up to five colours.

Stacking ten orbits of glue centres takes Sa from mean degree 9.94 to 18.32 at
6043 points, and every one of those carriers is 4-COLOURABLE in under a second
-- density alone never bought a chromatic number.  The tuned chain changes what
that density is worth: composing the forced pair with a rotated copy and tuning
the composite distance to exactly 1 makes the pair an EDGE, so the union of two
copies of the carrier refuses four without any spindle and keeps the carrier's
mean degree.

So this is a 5-chromatic unit-distance graph at degree 18.32, against a
previous ceiling of 13.82.  The per-vertex null model puts free@5 at about
1.2 % there, down from 16 % at degree 13.35, which is the regime where local
rigidity would have to start showing if it ever does.

With D^2 = 64/9 the tuning angle is cos phi = -119/128 and
sin phi = 3 sqrt247/128, already the carrier's own field, so none of this costs
a new radical.
"""
import sys, time, json, random, math
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 247))
rot60 = _rot60(F); Sa = build_Sa(F); g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
BEST = [265, 253, 211, 43, 139, 67, 241, 79, 277, 307]
cs = []
for b in BEST:
    for w in orbit(Sa[b]):
        if w not in cs: cs.append(w)
seen, U = set(Sa), list(Sa)
for w in cs:
    rot = g60.about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"carrier n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)

def cnf_for(gr, K, pin=True):
    N = gr.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(N)]
    for v in range(N):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in gr.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    if pin:
        for i, x in enumerate(gr.find_clique(3) or []):
            cl.append([X(x, i)])
            for c in range(K):
                if c != i:
                    cl.append([-X(x, c)])
    return cl, X

cl, X = cnf_for(g, 4, pin=False)
s = Solver(name="m22", bootstrap_with=cl)
assert s.solve()
rng = random.Random(23)
def read():
    p = set(l for l in s.get_model() if l > 0)
    return [next(c for c in range(4) if X(v, c) in p) for v in range(n)]
cols = [read()]
while len(cols) < 8:
    s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(4)])
    if not s.solve(): break
    cols.append(read())
buck = defaultdict(list)
for v in range(n):
    buck[tuple(c[v] for c in cols)].append(v)
DD = F.rational(Fr(64, 9)); pair = None
for vs in buck.values():
    for i, x in enumerate(vs):
        for y in vs[i+1:]:
            if (g.vertices[x] - g.vertices[y]).norm2() == DD and \
               not s.solve(assumptions=[X(x, 0), X(y, 1)]):
                pair = (x, y); break
        if pair: break
    if pair: break
s.delete()
assert pair, "no forced pair at 64/9"
v, q = pair
print(f"  forced pair {v},{q} confirmed   [{time.time()-t0:.0f}s]", flush=True)

COS = F.rational(Fr(9, 128) - 1)
SIN = F.sqrt(247) * F.rational(Fr(3, 128))
assert COS * COS + SIN * SIN == F.rational(1)
R = Rotation(COS, SIN)
V, Q = g.vertices[v], g.vertices[q]
def tau(p):
    d = Point(p.x - V.x, p.y - V.y); r = R(d)
    return Point(Q.x + r.x, Q.y + r.y)
assert (V - tau(Q)).norm2() == F.rational(1)
seen2, W = set(g.vertices), list(g.vertices)
for p in g.vertices:
    z = tau(p)
    if z not in seen2: seen2.add(z); W.append(z)
g2 = build_graph(W); n2 = g2.n; m2 = sum(len(a) for a in g2.adj) // 2
degs = [len(a) for a in g2.adj]
mu = sum(degs) / n2
pv = sum(4 * (0.75 ** d) for d in degs) / n2
print(f"  chained n={n2} m={m2} deg={mu:.2f}  per-vertex null free@5 "
      f"{100*pv:.3f} %   [{time.time()-t0:.0f}s]", flush=True)

cl4, _ = cnf_for(g2, 4)
t1 = time.time()
s4 = Solver(name="cd19", bootstrap_with=cl4); ok4 = s4.solve(); s4.delete()
print(f"  4-colourable: {ok4}   [{time.time()-t1:.0f}s]", flush=True)
assert not ok4
json.dump({"field_generators": list(F.gens), "n": n2, "m": m2,
           "mechanism": "ten stacked orbits, tuned chain to distance 1",
           "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                       [[c.numerator, c.denominator] for c in p.y.c]]
                      for p in g2.vertices]},
          open(f"{ROOT}/data/five_dense_10.json", "w"))
print("  written data/five_dense_10.json", flush=True)

cl5, X5 = cnf_for(g2, 5)
t1 = time.time()
s5 = Solver(name="cd19", bootstrap_with=cl5)
r5 = s5.solve()
print(f"  5-colourable: {r5}   [{time.time()-t1:.0f}s]", flush=True)
if not r5:
    print("  *** REFUSES FIVE ***", flush=True)
    json.dump({"field_generators": list(F.gens), "n": n2, "m": m2,
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g2.vertices]},
              open(f"{ROOT}/data/six_candidate.json", "w"))
else:
    pos = set(l for l in s5.get_model() if l > 0)
    col = [next(c for c in range(5) if X5(u, c) in pos) for u in range(n2)]
    free = sum(1 for u in range(n2)
               if len(set(col[w] for w in g2.adj[u])) < 4) / n2
    print(f"  free@5 = {100*free:.3f} %  against per-vertex null "
          f"{100*pv:.3f} %  ({free/max(pv,1e-12):.2f}x)", flush=True)
s5.delete()
