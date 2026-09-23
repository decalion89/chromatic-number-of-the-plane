"""Both tunings at once: a group and a chromatic number in the same object.

Tuning the DISTANCE moves a forced pair anywhere in [0, 2D].  Tuning the ANGLE
to 60 degrees makes the composition have order six, so the union of its copies
is C6-invariant and the whole orbit of the pivot is forced to one colour, on a
circle of radius R with R^2 = D^2.

Do both.  First move the forced pair to squared distance 1/3, which costs
cos phi1 = -125/128 and sin phi1 = sqrt(759)/128 with 759 = 3 * 11 * 23.  Then
tune the angle on THAT pair.  The orbit now lies on a circle of radius
1/sqrt3, and its six points are separated by 60 degrees, so the pairs two apart
are at distance

        2 R sin 60 = 2 * (1/sqrt3) * (sqrt3/2) = 1

exactly -- ADJACENT, and forced equal.  The union of the six copies therefore
has no 4-colouring at all, with no spindle anywhere in the construction, and it
carries C6 about a centre that belongs to none of the copies.

That is the object the whole line was aiming at: the group and the chromatic
number produced together, by tuning alone, from a single forced pair.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 23))
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
print(f"  carrier n={n}   [{time.time()-t0:.0f}s]", flush=True)

def cnf_for(gr, K):
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
    return cl, X

K = 4
cnf, X = cnf_for(g, K)
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
print(f"  forced pair at 64/9 confirmed   [{time.time()-t0:.0f}s]", flush=True)

# tuning one: move it to squared distance 1/3
C1 = F.rational(Fr(-125, 128)); S1 = F.sqrt(759) * F.rational(Fr(1, 128))
assert C1 * C1 + S1 * S1 == F.rational(1)
R1 = Rotation(C1, S1)
def tau1(p):
    d = Point(p.x - V.x, p.y - V.y); r = R1(d)
    return Point(Q.x + r.x, Q.y + r.y)
TQ = tau1(Q)
assert (V - TQ).norm2() == F.rational(Fr(1, 3)), "distance tuning missed"
seen2, W = set(g.vertices), list(g.vertices)
for p in g.vertices:
    z = tau1(p)
    if z not in seen2: seen2.add(z); W.append(z)
g2 = build_graph(W)
print(f"  after distance tuning n={g2.n}, pair at 1/3   "
      f"[{time.time()-t0:.0f}s]", flush=True)

# tuning two: angle 60 on the new pair
C2 = F.rational(Fr(1, 2)); S2 = F.sqrt(3) * F.rational(Fr(1, 2))
R2 = Rotation(C2, S2)
def tau2(p):
    d = Point(p.x - V.x, p.y - V.y); r = R2(d)
    return Point(TQ.x + r.x, TQ.y + r.y)
orbit, z = [], V
for _ in range(6):
    orbit.append(z); z = tau2(z)
assert z == V and len(set(orbit)) == 6, "tau2 must have order 6"
cnt = Counter(str((orbit[i] - orbit[j]).norm2())
              for i in range(6) for j in range(i + 1, 6))
print(f"  orbit squared distances: {dict(cnt)}", flush=True)
assert cnt.get("1", 0) > 0, "the orbit should contain unit distances"

seen3, Z = set(g2.vertices), list(g2.vertices)
cur = list(g2.vertices)
for _ in range(5):
    cur = [tau2(p) for p in cur]
    for p in cur:
        if p not in seen3: seen3.add(p); Z.append(p)
g3 = build_graph(Z); n3 = g3.n; m3 = sum(len(a) for a in g3.adj) // 2
pos = {p: i for i, p in enumerate(g3.vertices)}
inv = all(tau2(p) in pos for p in g3.vertices)
print(f"  union n={n3} m={m3} deg={2.0*m3/n3:.2f}  tau2-invariant {inv}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

cl3, _ = cnf_for(g3, 4)
t1 = time.time()
s3 = Solver(name="cd19", bootstrap_with=cl3); ok = s3.solve(); s3.delete()
print(f"  4-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
if not ok:
    print("  *** refuses four: a C6-invariant 5-chromatic graph, "
          "by tuning alone ***", flush=True)
    json.dump({"field_generators": list(F.gens), "n": n3, "m": m3,
               "mechanism": "distance tuned to 1/3, then angle tuned to order 6",
               "orbit": [pos[p] for p in orbit], "orbit_d2": dict(cnt),
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g3.vertices]},
              open(f"{ROOT}/data/five_twotune.json", "w"))
    print("  written data/five_twotune.json", flush=True)
