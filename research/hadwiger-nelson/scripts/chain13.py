"""Forced pairs compose, and the composite distance is tunable.

Every forced pair this project has ever produced sits at squared distance
64/9, 64/3, 256/9 or 16 -- distance 2.67 and up.  The disjunctive spindle needs
1/sqrt3 = 0.577.  The gadget is aimed at a regime where forcing does not occur,
and that, not the solver, is why every search returns zero.

But forcing is transitive and the distance is a free parameter.  Let H force
c(v) = c(q) with |v - q| = D, and let tau be the isometry taking v to q whose
rotational part is R_phi:

        tau(p) = q + R_phi (p - v)

Then tau(H) forces c(tau v) = c(tau q), i.e. c(q) = c(tau q), so H u tau(H)
forces c(v) = c(tau q) -- a NEW forced pair.  Its distance is

        v - tau(q) = (v - q) - R_phi (q - v) = u + R_phi u,   u = v - q

        |v - tau(q)|^2 = |u|^2 (2 + 2 cos phi) = 2 D^2 (1 + cos phi)

so the composite distance sweeps the whole interval [0, 2D] as phi turns, and
the angle for a target d is

        cos phi = d^2 / (2 D^2) - 1

which is RATIONAL whenever d^2 and D^2 are.  Only sin phi carries a radical,
and that radical is the entire cost of the move.

Here D^2 = 64/9 and d^2 = 1/3 give cos phi = -125/128 and
sin phi = sqrt(759)/128 with 759 = 3 * 11 * 23, so the composite pair lands at
exactly 1/sqrt3 in Q(sqrt3, sqrt5, sqrt7, sqrt11, sqrt23).

At that distance the ordinary spindle angle 2 arcsin(1/(2d)) is 120 degrees --
the same rotation the three-fold gadget uses -- so the pair can be spent at
once, and the result is a 5-chromatic graph in a field reached by tuning rather
than by luck.  The same recipe runs for any rational target, which turns the
field of a 5-chromatic unit-distance graph into a parameter.
"""
import sys, time, json, random
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
F = Field((3, 5, 7, 11, 23))
print(f"field {F.gens}, dim {F.dim}", flush=True)
rot60 = _rot60(F)
Sa = build_Sa(F)

def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for r in [Point(z.x, -z.y) for z in list(out)]:
        if r not in out: out.append(r)
    return out

seen, U = set(Sa), list(Sa)
for w in orb(Sa[25]):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        z = rot(p)
        if z not in seen: seen.add(z); U.append(z)
g = build_graph(U); n = g.n
print(f"carrier n={n} m={sum(len(a) for a in g.adj)//2}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            cnf.append([-X(v, a), -X(v, b)])
for x, y in g.edges():
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])
s = Solver(name="m22", bootstrap_with=cnf)
assert s.solve()
rng = random.Random(3)
def read():
    p = set(l for l in s.get_model() if l > 0)
    return [next(c for c in range(K) if X(v, c) in p) for v in range(n)]
cols = [read()]
while len(cols) < 10:
    s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    cols.append(read())
buck = defaultdict(list)
for v in range(n):
    buck[tuple(c[v] for c in cols)].append(v)
D2 = F.rational(Fr(64, 9))
pair = None
for vs in buck.values():
    for i, x in enumerate(vs):
        for y in vs[i+1:]:
            if (g.vertices[x] - g.vertices[y]).norm2() != D2:
                continue
            if not s.solve(assumptions=[X(x, 0), X(y, 1)]):
                pair = (x, y); break
        if pair: break
    if pair: break
s.delete()
assert pair, "no confirmed forced pair at 64/9"
v, q = pair
print(f"  forced pair {v},{q} at d^2 = 64/9, confirmed   "
      f"[{time.time()-t0:.0f}s]", flush=True)

# tune the composite distance to exactly 1/3
TARGET = Fr(1, 3)
cosphi = F.rational(TARGET / (2 * Fr(64, 9)) - 1)
sinphi = F.sqrt(759) * F.rational(Fr(1, 128))
assert cosphi * cosphi + sinphi * sinphi == F.rational(1)
print(f"  cos phi = {cosphi}, sin phi = sqrt759/128   (exact)", flush=True)
R = Rotation(cosphi, sinphi)
V, Q = g.vertices[v], g.vertices[q]
def tau(p):
    d = Point(p.x - V.x, p.y - V.y)
    r = R(d)
    return Point(Q.x + r.x, Q.y + r.y)
assert tau(V) == Q, "tau must carry v to q"
tq = tau(Q)
dd = (V - tq).norm2()
print(f"  |v - tau(q)|^2 = {dd}   (target 1/3)", flush=True)
assert dd == F.rational(TARGET)

seen2, W = set(g.vertices), list(g.vertices)
for p in g.vertices:
    z = tau(p)
    if z not in seen2: seen2.add(z); W.append(z)
g2 = build_graph(W); n2 = g2.n
pos = {p: i for i, p in enumerate(g2.vertices)}
print(f"  chained union n={n2} m={sum(len(a) for a in g2.adj)//2}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def cnf_for(gr, KK):
    N = gr.n
    Y = lambda u, c: 1 + u * KK + c
    cl = [[Y(u, c) for c in range(KK)] for u in range(N)]
    for u in range(N):
        for a in range(KK):
            for b in range(a + 1, KK):
                cl.append([-Y(u, a), -Y(u, b)])
    for x, y in gr.edges():
        for c in range(KK):
            cl.append([-Y(x, c), -Y(y, c)])
    return cl, Y

cl2, Y2 = cnf_for(g2, 4)
s2 = Solver(name="cd19", bootstrap_with=cl2)
iv, itq = pos[V], pos[tq]
forced = not s2.solve(assumptions=[Y2(iv, 0), Y2(itq, 1)])
print(f"  composite pair forced at four: {forced}   [{time.time()-t0:.0f}s]",
      flush=True)
s2.delete()
assert forced

# spend it: the spindle angle of d^2 = 1/3 is 120 degrees
rho = Rotation(F.rational(Fr(-1, 2)), F.sqrt(3) * F.rational(Fr(1, 2))).about(V)
assert (tq - rho(tq)).norm2() == F.rational(1)
seen3, Z = set(g2.vertices), list(g2.vertices)
for p in g2.vertices:
    z = rho(p)
    if z not in seen3: seen3.add(z); Z.append(z)
g3 = build_graph(Z); n3 = g3.n; m3 = sum(len(a) for a in g3.adj) // 2
print(f"  spindled n={n3} m={m3} deg={2.0*m3/n3:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
cl3, Y3 = cnf_for(g3, 4)
t1 = time.time()
s3 = Solver(name="cd19", bootstrap_with=cl3); ok = s3.solve(); s3.delete()
print(f"  4-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
if not ok:
    print("  *** refuses four: 5-chromatic, by tuning ***", flush=True)
    json.dump({"field_generators": list(F.gens), "n": n3, "m": m3,
               "mechanism": "composed forcing tuned to d^2 = 1/3, then spindled",
               "cos_phi": [-125, 128], "sin_phi": "sqrt(759)/128",
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g3.vertices]},
              open(f"{ROOT}/data/five_tuned_third.json", "w"))
    print("  written data/five_tuned_third.json", flush=True)
