"""Spend the free angle on the orbit the angle-tuning manufactured.

The order-six tuning turns any forced pair at 64/9 into a C6-invariant carrier
carrying a forced ORBIT of six points, at squared distances 64/9, 64/3 and
256/9 from one another.  From the pivot orbit[0] the others sit at

    orbit[1], orbit[5]  ->  64/9        orbit[2], orbit[4]  ->  64/3
    orbit[3]            -> 256/9

and 64/3 with 256/9 is the one pair whose circles meet, |sqrt r1 - sqrt r2| =
0.715 <= 1, so it is the one pair admitting a rotation whose angle is not
forced to make 4d^2 - 1 a square.

With a = orbit[0], b = orbit[2], c = orbit[3], u = b - a and v = c - a, any
rotation carrying b to distance exactly 1 from c closes the argument, because
a, b, c share a colour and rho fixes a:

    |rho(u) - v|^2 = r1 + r2 - 2 <rho(u), v>,   <rho(u), v> = cos t P + sin t Q
    P = u.v,  Q = u x v,  R = (r1 + r2 - 1)/2,  P^2 + Q^2 = r1 r2

one LINEAR equation beside cos^2 + sin^2 = 1, so both roots are explicit:

    cos t = (P R +- Q D)/(r1 r2),   sin t = (Q R -+ P D)/(r1 r2)
    D = sqrt(r1 r2 - R^2) = sqrt(16384/27 - (439/18)^2) = sqrt(3887)/18
                          = 13 sqrt(23) / 18

so the whole construction costs one radical, sqrt23, on top of the carrier's
own field.  The result is a 5-chromatic unit-distance graph in
Q(sqrt3, sqrt11, sqrt23, sqrt247) reached entirely by the new recipe: tune the
angle to order six, then spend the free angle the orbit provides.
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
F = Field((3, 11, 23, 247))
print(f"field {F.gens}, dim {F.dim}", flush=True)
rot60 = _rot60(F); Sa = build_Sa(F)
def orb6(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for r in [Point(z.x, -z.y) for z in list(out)]:
        if r not in out: out.append(r)
    return out
seen, U = set(Sa), list(Sa)
for w in orb6(Sa[25]):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        z = rot(p)
        if z not in seen: seen.add(z); U.append(z)
g = build_graph(U); n = g.n
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
s = Solver(name="m22", bootstrap_with=cnf); assert s.solve()
rng = random.Random(6)
def read():
    p = set(l for l in s.get_model() if l > 0)
    return [next(c for c in range(K) if X(v, c) in p) for v in range(n)]
cols = [read()]
while len(cols) < 8:
    s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
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
assert pair
V, Q = g.vertices[pair[0]], g.vertices[pair[1]]
COS = F.rational(Fr(1, 2)); SIN = F.sqrt(3) * F.rational(Fr(1, 2))
RR = Rotation(COS, SIN)
def tau(p):
    d = Point(p.x - V.x, p.y - V.y); r = RR(d)
    return Point(Q.x + r.x, Q.y + r.y)
orbit, z = [], V
for _ in range(6):
    orbit.append(z); z = tau(z)
assert z == V and len(set(orbit)) == 6
seen2, W = set(g.vertices), list(g.vertices)
cur = list(g.vertices)
for _ in range(5):
    cur = [tau(p) for p in cur]
    for p in cur:
        if p not in seen2: seen2.add(p); W.append(p)
g2 = build_graph(W); pos = {p: i for i, p in enumerate(g2.vertices)}
print(f"  C6 carrier n={g2.n} m={sum(len(a) for a in g2.adj)//2}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

A_, B_, C_ = orbit[0], orbit[2], orbit[3]
assert (A_ - B_).norm2() == F.rational(Fr(64, 3))
assert (A_ - C_).norm2() == F.rational(Fr(256, 9))
u = Point(B_.x - A_.x, B_.y - A_.y)
vv = Point(C_.x - A_.x, C_.y - A_.y)
P = u.x * vv.x + u.y * vv.y
Qq = u.x * vv.y - u.y * vv.x
Rr = F.rational((Fr(64, 3) + Fr(256, 9) - 1) / 2)
Dd = F.sqrt(23) * F.rational(Fr(13, 18))
inv = F.rational(Fr(27, 16384))
built = []
for sgn in (1, -1):
    co = (P * Rr + Qq * Dd * F.rational(sgn)) * inv
    si = (Qq * Rr - P * Dd * F.rational(sgn)) * inv
    if co * co + si * si == F.rational(1):
        r = Rotation(co, si)
        if (r.about(A_)(B_) - C_).norm2() == F.rational(1):
            built.append(r)
print(f"  {len(built)} free-angle rotations, exact   [{time.time()-t0:.0f}s]",
      flush=True)
assert built

rho = built[0].about(A_)
seen3, Z = set(g2.vertices), list(g2.vertices)
for p in g2.vertices:
    z = rho(p)
    if z not in seen3: seen3.add(z); Z.append(z)
g3 = build_graph(Z); n3 = g3.n; m3 = sum(len(a) for a in g3.adj) // 2
print(f"  spindled n={n3} m={m3} deg={2.0*m3/n3:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
Y = lambda u_, c: 1 + u_ * K + c
cl = [[Y(u_, c) for c in range(K)] for u_ in range(n3)]
for u_ in range(n3):
    for a in range(K):
        for b in range(a + 1, K):
            cl.append([-Y(u_, a), -Y(u_, b)])
for x, y in g3.edges():
    for c in range(K):
        cl.append([-Y(x, c), -Y(y, c)])
t1 = time.time()
s3 = Solver(name="cd19", bootstrap_with=cl); ok = s3.solve(); s3.delete()
print(f"  4-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
if not ok:
    print("  *** refuses four: 5-chromatic by the full new recipe ***",
          flush=True)
    json.dump({"field_generators": list(F.gens), "n": n3, "m": m3,
               "mechanism": "angle tuned to order 6, then the free angle spent",
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g3.vertices]},
              open(f"{ROOT}/data/five_order6_free.json", "w"))
    print("  written data/five_order6_free.json", flush=True)
