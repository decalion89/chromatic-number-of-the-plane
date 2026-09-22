"""Minimise a symmetric carrier without breaking its symmetry.

The spindle of the 1021-point symmetric carrier has 7141 vertices, and every
5-colouring of it costs cadical minutes, so the forced-pair test at five --
the only question that matters -- would take a day.  Shrinking the carrier
shrinks the spindle proportionally, but ordinary minimisation deletes one
vertex at a time and destroys the group, and the group is what makes the test
affordable: it multiplies the colourings twelve-fold and divides the
confirmations by twelve.

So delete ORBITS instead of vertices.  An orbit is at most twelve points, the
carrier stays invariant after every deletion, and the pair being kept forced is
itself an orbit, so it survives as a whole.  The result is a smaller carrier
with the same group and the same forcing.
"""
import sys, time, random, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
VI = int(sys.argv[1]) if len(sys.argv) > 1 else 25
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
def orb(p, refl=True):
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
for w in orb(Sa[VI]):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
pos = {p: i for i, p in enumerate(g.vertices)}
# orbits of the point set under the rotation group (the reflection may not be
# a symmetry of the spindled object later, so keep to the rotations)
orbits, mark = [], [False] * n
for v in range(n):
    if mark[v]: continue
    o, q = [], g.vertices[v]
    for _ in range(6):
        i = pos[q]
        if not mark[i]: mark[i] = True; o.append(i)
        q = rot60(q)
    orbits.append(o)
print(f"carrier n={n}, {len(orbits)} rotation orbits   [{time.time()-t0:.0f}s]",
      flush=True)
K = 4
X = lambda v, c: 1 + v * K + c
SEL = lambda v: 1 + n * K + v
base = [[-SEL(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        base.append([-X(u, c), -X(v, c)])
# find one forced pair to keep
cnf0 = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf0.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf0); s.solve()
p0 = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in p0) for v in range(n)]]
rng = random.Random(SEED)
while len(cols) < 8:
    s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
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
            s2 = Solver(name="m22", bootstrap_with=cnf0)
            d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
            if not d:
                dd = (U[a] - U[b]).norm2()
                if dd.is_rational():
                    pair = (a, b, Fr(dd.c[0])); break
        if pair: break
    if pair: break
A_, B_, D2 = pair
print(f"  keeping the forcing of v{A_},v{B_} at d^2={D2}   [{time.time()-t0:.0f}s]",
      flush=True)
keepset = set()
for o in orbits:
    if A_ in o or B_ in o: keepset.update(o)
cnf = list(base) + [[X(A_, 0)], [X(B_, 1)]]
def forces(sub):
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[SEL(v) for v in sorted(sub)]); s.delete()
    return not ok
cur = set(range(n))
assert forces(cur)
order = list(range(len(orbits)))
rng.shuffle(order)
dropped = 0
for oi in order:
    o = orbits[oi]
    if any(v in keepset for v in o): continue
    trial = cur - set(o)
    if forces(trial):
        cur = trial; dropped += 1
print(f"  dropped {dropped} orbits -> carrier of {len(cur)} points"
      f"   [{time.time()-t0:.0f}s]", flush=True)
V = [U[v] for v in sorted(cur)]
gv = build_graph(V)
S = set(gv.vertices)
print(f"  still rot60-invariant: {all(rot60(p) in S for p in gv.vertices)}",
      flush=True)
rot = rotation_joining(D2, F)
piv = orb(U[A_], refl=False)
seen2, W = set(V), list(V)
for w in piv:
    r = rot.about(w)
    for p in V:
        z = r(p)
        if z not in seen2: seen2.add(z); W.append(z)
gw = build_graph(W); mw = sum(len(a) for a in gw.adj) // 2
print(f"  spindled over {len(piv)} pivots: n={gw.n} m={mw} deg={2.0*mw/gw.n:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
tri = gw.find_clique(3)
for KK in (4, 5):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(gw.n)]
    for u, v in gw.edges():
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
    if KK == 5: sys.exit()
json.dump({"field_generators": list(F.gens), "n": gw.n, "m": mw,
           "carrier": len(cur), "d2": str(D2),
           "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                       [[c.numerator, c.denominator] for c in p.y.c]]
                      for p in gw.vertices]},
          open(f"/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/equimin_{VI}_{SEED}.json", "w"))
print("  saved", flush=True)
