"""The other forced orbit, spindled alone: a third field.

The symmetric carrier's 153 forced pairs fall in two orbits.  The first, at
d^2 = 64/9, has been spent: it gives the 7141-vertex C6-invariant graph in
Q(sqrt3,sqrt11,sqrt247).  The second, nine pairs at d^2 = 64/3, needs a
different radical --

    cos t = 1 - (1/2)/(64/3) = 125/128,  sin t = sqrt(759)/128
    759 = 3 * 11 * 23

so sqrt(759) = sqrt3 sqrt11 sqrt23 and the graph lives in Q(sqrt3,sqrt11,
sqrt23), a field this project has never produced a graph in.  Spending it
alone, over its whole orbit of pivots, gives a different 5-chromatic graph
rather than a bigger one.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 23))
print(f"field {F.gens}, dimension {F.dim}", flush=True)
rot60 = _rot60(F); Sa = build_Sa(F)
def orb(p, refl=False):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if refl:
        for q in list(out):
            r = Point(q.x, -q.y)
            if r not in out: out.append(r)
    return out
t0 = time.time()
seen, U = set(Sa), list(Sa)
for w in orb(Sa[25], True):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
print(f"carrier n={n} m={sum(len(a) for a in g.adj)//2}   [{time.time()-t0:.0f}s]",
      flush=True)
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf); s.solve()
p0 = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in p0) for v in range(n)]]
rng = random.Random(7)
while len(cols) < 16:
    s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 40)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in cols)].append(v)
pair = None
for vs in buck.values():
    for i, a in enumerate(vs):
        for b in vs[i+1:]:
            s2 = Solver(name="m22", bootstrap_with=cnf)
            d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
            if not d:
                dd = (U[a] - U[b]).norm2()
                if dd.is_rational() and Fr(dd.c[0]) == Fr(64, 3):
                    pair = (a, b); break
        if pair: break
    if pair: break
if pair is None:
    print("  no forced pair at 64/3 found in this sample", flush=True); sys.exit()
a, b = pair
print(f"  forced pair v{a},v{b} at d^2=64/3   [{time.time()-t0:.0f}s]", flush=True)
spin = rotation_joining(Fr(64, 3), F)
seen2, V = set(U), list(U)
for w in orb(U[a]):
    rot = spin.about(w)
    for p in U:
        z = rot(p)
        if z not in seen2: seen2.add(z); V.append(z)
gv = build_graph(V); nv = gv.n; mv = sum(len(x) for x in gv.adj) // 2
print(f"  spindled over 6 pivots: n={nv} m={mv} deg={2.0*mv/nv:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
tri = gv.find_clique(3)
for KK in (4, 5):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(nv)]
    for u, v in gv.edges():
        for c in range(KK):
            cc.append([-XX(u, c), -XX(v, c)])
    for i, v in enumerate(tri):
        cc.append([XX(v, i)])
        for c in range(KK):
            if c != i: cc.append([-XX(v, c)])
    t1 = time.time()
    sv = Solver(name="cd19", bootstrap_with=cc); ok = sv.solve(); sv.delete()
    print(f"  {KK}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if ok: break
    print(f"  *** refuses {KK} ***", flush=True)
    if KK == 4:
        json.dump({"field_generators": list(F.gens), "n": nv, "m": mv,
                   "recipe": ["Sa glued at the D6 orbit of Sa[25]",
                              "spindled at d^2 = 64/3 over the C6 orbit of a pivot",
                              "cos 125/128, sin sqrt(759)/128, 759 = 3*11*23"],
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gv.vertices]},
                  open("/home/user/darwin-50/research/hadwiger-nelson/data/five_23.json", "w"))
        print("  written data/five_23.json", flush=True)
