"""Tune the angle instead of the distance, and the chain closes into a group.

The composition tau(p) = q + R_phi (p - v) is a rotation -- its linear part is
one -- so it is a rotation by phi about the point w fixed by it, and since
tau(v) = q the two points v and q lie on a circle about w separated by phi:

        |v - q|^2 = 2 R^2 (1 - cos phi),      R = |w - v|

Earlier the angle was chosen to place the COMPOSITE DISTANCE.  Choose it to fix
the ORDER instead.  phi = 60 degrees needs only sqrt3, and then tau^6 = id, so

  * the union of the six copies tau^j(H) is C6-INVARIANT about w -- a group
    the carrier did not have, about a centre that is not in it;
  * the forcing chains all the way round.  In tau^j(H) the pair (tau^j v,
    tau^j q) is forced, and tau^j v = tau^{j-1} q, so every consecutive pair of
    the ORBIT of v is forced equal and the whole orbit takes one colour;
  * R^2 = D^2/(2(1 - cos 60)) = D^2, so with D^2 = 64/9 the orbit is six points
    on a circle of radius 8/3, and their mutual squared distances are
    64/9, 64/3 and 256/9 -- three different radii, forced, in one class.

That last line is the point.  The "class of three at different radii" that
turned up by luck on one carrier, and that the free-angle spindle needs, is
what this tuning produces from ANY forced pair at 64/9.  64/3 with 256/9 gives
|sqrt(r1) - sqrt(r2)| = 0.71 <= 1, so the circles meet and the free angle
exists; 64/9 with either gives 1.95 and nothing.
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
F = Field((3, 11, 247))
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
assert pair
v, q = pair
V, Q = g.vertices[v], g.vertices[q]
print(f"  forced pair at 64/9 confirmed   [{time.time()-t0:.0f}s]", flush=True)

# phi = 60 degrees
COS = F.rational(Fr(1, 2)); SIN = F.sqrt(3) * F.rational(Fr(1, 2))
R = Rotation(COS, SIN)
def tau(p):
    d = Point(p.x - V.x, p.y - V.y); r = R(d)
    return Point(Q.x + r.x, Q.y + r.y)
assert tau(V) == Q
z = V
orbit = []
for _ in range(6):
    orbit.append(z); z = tau(z)
assert z == V, "tau must have order 6 on v"
assert len(set(orbit)) == 6
print(f"  tau has order 6; orbit of v is {len(set(orbit))} points", flush=True)
ds = sorted({str((orbit[i] - orbit[j]).norm2())
             for i in range(6) for j in range(i + 1, 6)})
print(f"  orbit squared distances: {ds}", flush=True)

seen2, W = set(g.vertices), list(g.vertices)
cur = list(g.vertices)
for _ in range(5):
    cur = [tau(p) for p in cur]
    for p in cur:
        if p not in seen2: seen2.add(p); W.append(p)
g2 = build_graph(W); n2 = g2.n; m2 = sum(len(a) for a in g2.adj) // 2
pos = {p: i for i, p in enumerate(g2.vertices)}
inv = all(tau(p) in pos for p in g2.vertices)
print(f"  union n={n2} m={m2} deg={2.0*m2/n2:.2f}  tau-invariant: {inv}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

Y = lambda u, c: 1 + u * K + c
cl = [[Y(u, c) for c in range(K)] for u in range(n2)]
for u in range(n2):
    for a in range(K):
        for b in range(a + 1, K):
            cl.append([-Y(u, a), -Y(u, b)])
for x, y in g2.edges():
    for c in range(K):
        cl.append([-Y(x, c), -Y(y, c)])
s2 = Solver(name="cd19", bootstrap_with=cl)
ok = s2.solve()
print(f"  4-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if ok:
    idx = [pos[p] for p in orbit]
    bad = []
    for i in range(6):
        for j in range(i + 1, 6):
            if s2.solve(assumptions=[Y(idx[i], 0), Y(idx[j], 1)]):
                bad.append((i, j))
    print(f"  orbit forced-equal: {'ALL 15 PAIRS' if not bad else bad}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"field_generators": list(F.gens), "n": n2, "m": m2,
               "mechanism": "tuned to angle 60, tau of order 6",
               "forced_orbit": idx, "orbit_d2": ds,
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g2.vertices]},
              open(f"{ROOT}/data/order6_carrier.json", "w"))
    print("  written data/order6_carrier.json", flush=True)
else:
    print("  *** the union already refuses four ***", flush=True)
s2.delete(); s.delete()
