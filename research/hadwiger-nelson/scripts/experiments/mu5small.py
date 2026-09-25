"""The same decisive test, on the smallest dense 5-chromatic object available.

five_dense_2 is 6925 points and its first mu call has been running half an
hour, so the same question is asked of the 1-orbit carrier chained to distance
1: 5378 points at mean degree 15.92, the exact five-colour analogue of the
2689-point carrier whose mu_4 saturated at 3 on all 150 of its richest points.

Only "is mu above 2?" is asked, which is a single call per candidate: forbid
three colours on N(p) and see whether a colouring survives.  A yes is mu <= 2
and the candidate is done; a NO is the first time anything in this project has
exceeded 2 at five colours.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
F = Field((3, 11, 247))
one = F.rational(Fr(1)); half = F.rational(Fr(1, 2))
rot60 = _rot60(F); r30 = Rotation(F.sqrt(3) * half, half)
Sa = build_Sa(F); g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
seen, U = set(Sa), list(Sa)
for w in orbit(Sa[265]):
    rot = g60.about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
print(f"carrier n={n} deg={2.0*sum(len(a) for a in g.adj)/2/n:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
# forced pair at 64/9, then the chain tuned to distance 1
K4 = 4
X4 = lambda v, c: 1 + v * K4 + c
cnf = [[X4(v, c) for c in range(K4)] for v in range(n)]
for v in range(n):
    for a in range(K4):
        for b in range(a + 1, K4):
            cnf.append([-X4(v, a), -X4(v, b)])
for x, y in g.edges():
    for c in range(K4):
        cnf.append([-X4(x, c), -X4(y, c)])
s = Solver(name="m22", bootstrap_with=cnf); assert s.solve()
import random
rng = random.Random(9)
def read():
    p = set(l for l in s.get_model() if l > 0)
    return [next(c for c in range(K4) if X4(v, c) in p) for v in range(n)]
cols = [read()]
while len(cols) < 8:
    s.add_clause([-X4(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X4(v, c)
                  for v in range(n) for c in range(K4)])
    if not s.solve(): break
    cols.append(read())
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in cols)].append(v)
DD = F.rational(Fr(64, 9)); pair = None
for vs in buck.values():
    for i, x in enumerate(vs):
        for y in vs[i+1:]:
            if (g.vertices[x] - g.vertices[y]).norm2() == DD and \
               not s.solve(assumptions=[X4(x, 0), X4(y, 1)]):
                pair = (x, y); break
        if pair: break
    if pair: break
s.delete()
assert pair
V, Q = g.vertices[pair[0]], g.vertices[pair[1]]
COS = F.rational(Fr(9, 128) - 1); SIN = F.sqrt(247) * F.rational(Fr(3, 128))
assert COS * COS + SIN * SIN == one
R = Rotation(COS, SIN)
def tau(p):
    d = Point(p.x - V.x, p.y - V.y); r = R(d)
    return Point(Q.x + r.x, Q.y + r.y)
assert (V - tau(Q)).norm2() == one
seen2, W = set(g.vertices), list(g.vertices)
for p in g.vertices:
    z = tau(p)
    if z not in seen2: seen2.add(z); W.append(z)
gg = build_graph(W); N = gg.n
print(f"chained n={N} deg={2.0*sum(len(a) for a in gg.adj)/2/N:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
S = set(gg.vertices)
cells = defaultdict(list)
for i, p in enumerate(gg.vertices):
    cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
gx = [q.fx for q in gg.vertices]; gy = [q.fy for q in gg.vertices]
deg_order = sorted(range(N), key=lambda v: -len(gg.adj[v]))
cand = set()
for c in deg_order[:20]:
    for rot in (rot60.about(gg.vertices[c]), r30.about(gg.vertices[c])):
        for p in gg.vertices:
            z = rot(p)
            if z not in S: cand.add(z)
scored, seenn = [], set()
for z in cand:
    zx, zy = z.fx, z.fy
    cx, cy = int(zx // 1), int(zy // 1)
    nb = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for i in cells.get((cx + dx, cy + dy), ()):
                ex, ey = gx[i] - zx, gy[i] - zy
                if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                   (gg.vertices[i] - z).norm2() == one:
                    nb.append(i)
    if len(nb) >= 10:
        t = tuple(sorted(nb))
        if t not in seenn:
            seenn.add(t); scored.append((len(nb), t))
scored.sort(key=lambda u: -u[0])
print(f"  {len(scored)} neighbourhoods >= 10, largest "
      f"{scored[0][0] if scored else 0}   [{time.time()-t0:.0f}s]", flush=True)
hist = Counter()
with MuSolver(gg, k=5, budget=3_000_000) as ms:
    print(f"  base colouring done, colourable={ms.colourable}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    for sz, nb in scored[:40]:
        r = ms.at_most_two(list(nb))
        hist[{True: 2, False: 3, None: -1}[r]] += 1
        if r is False:
            v = ms.mu(list(nb))
            print(f"  *** mu_5 = {v} with |N|={sz} -- ABOVE TWO AT FIVE ***",
                  flush=True)
        if sum(hist.values()) % 5 == 0:
            print(f"    ..{sum(hist.values())}: {dict(hist)}"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"  mu_5 over {sum(hist.values())}: {dict(hist)}   "
      f"(2 = at most two, 3 = above, -1 = budget)   [{time.time()-t0:.0f}s]",
      flush=True)
