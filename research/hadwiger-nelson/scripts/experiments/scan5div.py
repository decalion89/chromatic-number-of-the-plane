"""Forced-SAME pairs at the distances a five-denominator rotation can close.

In an integral multiquadratic graph every coset colouring survives any
integral in-field rotation, so a forced pair whose closing rotation is
integral at 5 cannot exist (its spindle would be 6-chromatic and coset
5-colourable at once).  The distances whose closure carries 5 in the
denominator, and lies in the Moser field, are
    d^2 = 5/9   closed by (1 + 3 sqrt-11)/10     cos 1/10
    d^2 = 25    closed by (49 + 3 sqrt-11)/50    cos 49/50  (Exoo-Ismailescu)
Scan exactly those classes: sampled colourings discard, the solver settles.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict, deque
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver
ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1]; SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 40
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n; E = list(G.edges())
adj = [[] for _ in range(n)]
for a, b in E: adj[a].append(b); adj[b].append(a)
hx = np.array([q.fx for q in G.vertices]); hy = np.array([q.fy for q in G.vertices])
TARGETS = {Fr(5, 9): "cos 1/10", Fr(25): "cos 49/50"}
pairs = []
cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 5.5), int(hy[i] // 5.5))].append(i)
for i in range(n):
    cx, cy = int(hx[i] // 5.5), int(hy[i] // 5.5)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i: continue
                dd = (hx[i] - hx[j])**2 + (hy[i] - hy[j])**2
                for T in TARGETS:
                    if abs(dd - float(T)) < 1e-9 and G.vertices[i].dist2(G.vertices[j]) == F.rational(T):
                        pairs.append((i, j, T))
print(f"{NAME}: n={n} m={len(E)}; pairs: " +
      ", ".join(f"d^2={T}: {sum(1 for p in pairs if p[2] == T)}" for T in TARGETS) +
      f"   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
s = Solver(name="cd19")
for v in range(n): s.add_clause([X(v, c) for c in range(K)])
for a, b in E:
    for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
def col_of():
    m = s.get_model()
    return np.array([next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)])
random.seed(3)
def kempe(col):
    v = random.randrange(n); a = col[v]; b = random.choice([c for c in range(K) if c != a])
    comp = {v}; q = deque([v])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in comp and (col[w] == a or col[w] == b): comp.add(w); q.append(w)
    idx = np.fromiter(comp, int); ca = col[idx] == a; col[idx[ca]] = b; col[idx[~ca]] = a
alive = list(range(len(pairs)))
PI = np.array([p[0] for p in pairs]); PJ = np.array([p[1] for p in pairs])
seen_same = np.zeros(len(pairs), int); nsamp = 0
for t in range(SAMPLES):
    v = random.randrange(n); c = random.randrange(K)
    assert s.solve(assumptions=[X(v, c)]); col = col_of()
    for _ in range(200):
        kempe(col)
        alive = [k for k in alive if col[PI[k]] == col[PJ[k]]]
        seen_same += (col[PI] == col[PJ]); nsamp += 1
    if t % 10 == 9:
        print(f"    {t+1} solver samples: {len(alive)} pairs never split   [{time.time()-t0:.0f}s]", flush=True)
for T in TARGETS:
    ks = [k for k in range(len(pairs)) if pairs[k][2] == T]
    if ks:
        ps = seen_same[ks] / nsamp
        print(f"  d^2={T}: mean P(same)={ps.mean():.3f}, max {ps.max():.3f}", flush=True)
forced = []; calls = 0
nsel = n * K + 1
while alive:
    k = alive[0]; i, j, T = pairs[k]
    sel = nsel; nsel += 1
    for c in range(K): s.add_clause([-sel, -X(i, c), -X(j, c)])
    r = s.solve(assumptions=[sel]); calls += 1
    c2 = col_of() if r else None
    s.add_clause([-sel])
    if not r:
        forced.append((i, j, str(T)))
        print(f"  *** FORCED SAME AT d^2={T} ({TARGETS[T]} closes it): ({i},{j}) -- VERIFY ***", flush=True)
        json.dump({"graph": NAME, "forced_same": forced}, open(f"{ROOT}/data/forced5div_{NAME}", "w"))
        alive = alive[1:]; continue
    for _ in range(200): kempe(c2)
    alive = [q for q in alive if c2[PI[q]] == c2[PJ[q]]]
    print(f"    solver call {calls}: {len(alive)} left   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  forced-same pairs at five-closable distances: {len(forced)}   [{time.time()-t0:.0f}s]", flush=True)
