"""Forced-DIFFERENT pairs at distance 2 -- the gadget Exoo-Ismailescu makes decisive.

Every 5-colouring of the plane has a monochromatic pair at distance 1 or 2
(Exoo-Ismailescu 2020).  So one unit-distance graph H with two points u, v a
distance 2 apart that NO 5-colouring of H colours alike proves chi(R^2) >= 6:
colour the plane with five colours, every congruent copy of H splits every
2-pair, and the colouring would then be a proper colouring of the {1,2}-graph.

The scan: every exact distance-2 pair of the graph is a candidate; any colouring
that gives it one colour discards it.  Kempe walks supply the colourings, and
every survivor gets the real test -- one call assuming c(u) = c(v).

Alongside, the same walk measures how often each distance class is monochromatic,
so that if distance 2 is not the easiest place to force a split, the data says
which distance is.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict, deque
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "tight_hexagon_4159.json"
SWAPS = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n
E = list(G.edges())
adj = [[] for _ in range(n)]
for a, b in E: adj[a].append(b); adj[b].append(a)
hx = np.array([q.fx for q in G.vertices]); hy = np.array([q.fy for q in G.vertices])
print(f"{NAME}: n={n} edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            cnf.append([-X(v, a), -X(v, b)])
for a, b in E:
    for c in range(K):
        cnf.append([-X(a, c), -X(b, c)])
s = Solver(name="cd19", bootstrap_with=cnf)
assert s.solve(), "not 5-colourable"
def colouring():
    m = s.get_model()
    return np.array([next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)])
col = colouring()
print(f"  first colouring   [{time.time()-t0:.0f}s]", flush=True)

# all pairs closer than 3, by squared distance class
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 3), int(hy[i] // 3))].append(i)
I, J = [], []
for i in range(n):
    cx, cy = int(hx[i] // 3), int(hy[i] // 3)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j > i: I.append(i); J.append(j)
I = np.array(I); J = np.array(J)
D2 = (hx[I] - hx[J])**2 + (hy[I] - hy[J])**2
keep = (D2 < 9.0) & (np.abs(D2 - 1.0) > 1e-9)
I, J, D2 = I[keep], J[keep], D2[keep]
cls = np.round(D2, 7)
print(f"  {len(I)} non-edge pairs closer than 3   [{time.time()-t0:.0f}s]", flush=True)
two = np.where(np.abs(D2 - 4.0) < 1e-9)[0]
four = F.rational(4)
two = np.array([t for t in two if G.vertices[I[t]].dist2(G.vertices[J[t]]) == four])
print(f"  {len(two)} exact distance-2 pairs   [{time.time()-t0:.0f}s]", flush=True)
alive = np.ones(len(two), bool)

random.seed(1)
def kempe(col):
    v = random.randrange(n)
    a = col[v]; b = random.choice([c for c in range(K) if c != a])
    comp = {v}; q = deque([v])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in comp and (col[w] == a or col[w] == b):
                comp.add(w); q.append(w)
    idx = np.fromiter(comp, int)
    ca = col[idx] == a
    col[idx[ca]] = b; col[idx[~ca]] = a
same = np.zeros(len(I)); samples = 0
def absorb(col):
    global samples, same
    eq = col[I] == col[J]
    same += eq; samples += 1
    alive[eq[two]] = False
absorb(col)
for t in range(1, SWAPS + 1):
    kempe(col)
    if t % 20 == 0: absorb(col)
    if t % 5000 == 0:
        print(f"    {t} swaps: {alive.sum()} distance-2 pairs never yet monochromatic"
              f"   [{time.time()-t0:.0f}s]", flush=True)
assert all(col[a] != col[b] for a, b in E)
# distance statistics
by = defaultdict(lambda: [0.0, 0])
for c, sm in zip(cls, same / samples):
    by[c][0] += sm; by[c][1] += 1
rows = sorted(((v[0] / v[1], c, v[1]) for c, v in by.items() if v[1] >= 200))
print("  least often monochromatic distance classes (d^2, pairs, P(same)):")
for p, c, m in rows[:15]: print(f"    d^2={c:.6f}  pairs={m:6d}  P(same)={p:.4f}")
p2 = by[4.0]; print(f"  distance 2: pairs={p2[1]} P(same)={p2[0]/max(p2[1],1):.4f}")
print(f"  (a random pair would be 0.2)")

cand = [int(t) for t in two[alive]]
print(f"  survivors after the walk: {len(cand)}   [{time.time()-t0:.0f}s]", flush=True)
nsel = n * K + 1; forced = []; calls = 0
while cand:
    t = cand[0]; i, j = int(I[t]), int(J[t])
    sel = nsel; nsel += 1
    for c in range(K):
        s.add_clause([-sel, -X(i, c), X(j, c)]); s.add_clause([-sel, X(i, c), -X(j, c)])
    r = s.solve(assumptions=[sel]); calls += 1
    c2 = colouring() if r else None
    s.add_clause([-sel])
    if not r:
        forced.append((i, j))
        print(f"  *** FORCED APART AT FIVE, DISTANCE 2: ({i},{j}) -- VERIFY ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        json.dump({"graph": NAME, "forced_apart_distance_2": forced},
                  open(f"{ROOT}/data/dist2_forced.json", "w"))
        cand = cand[1:]; continue
    for _ in range(300): kempe(c2)
    cand = [t for t in cand if c2[I[t]] != c2[J[t]]]
    print(f"    solver call {calls}: {len(cand)} left   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  distance-2 pairs forced apart at five: {len(forced)}   [{time.time()-t0:.0f}s]")
