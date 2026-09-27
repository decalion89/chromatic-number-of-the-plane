"""P(same colour) by distance class over Kempe-sampled 5-colourings, full table."""
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
NAME = sys.argv[1]; SWAPS = int(sys.argv[2]); K = 5
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G = build_graph(P); n = G.n; E = list(G.edges())
adj = [[] for _ in range(n)]
for a, b in E: adj[a].append(b); adj[b].append(a)
hx = np.array([q.fx for q in G.vertices]); hy = np.array([q.fy for q in G.vertices])
X = lambda v, c: 1 + v * K + c
s = Solver(name="cd19")
for v in range(n): s.add_clause([X(v, c) for c in range(K)])
for a, b in E:
    for c in range(K): s.add_clause([-X(a, c), -X(b, c)])
assert s.solve(); m = s.get_model()
col = np.array([next(c for c in range(K) if m[X(v, c) - 1] > 0) for v in range(n)])
I, J = np.triu_indices(n, 1)
D2 = (hx[I] - hx[J])**2 + (hy[I] - hy[J])**2
k = D2 < 9.0; I, J, D2 = I[k], J[k], D2[k]; cls = np.round(D2, 6)
random.seed(2)
def kempe(col):
    v = random.randrange(n); a = col[v]; b = random.choice([c for c in range(K) if c != a])
    comp = {v}; q = deque([v])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in comp and (col[w] == a or col[w] == b): comp.add(w); q.append(w)
    idx = np.fromiter(comp, int); ca = col[idx] == a; col[idx[ca]] = b; col[idx[~ca]] = a
same = np.zeros(len(I)); S = 0
for t in range(SWAPS):
    kempe(col)
    if t % 20 == 0: same += (col[I] == col[J]); S += 1
u, inv = np.unique(cls, return_inverse=True)
tot = np.bincount(inv, weights=same / S); cnt = np.bincount(inv)
want = [5/9, 25, 5/3, 5, 20/9, 25/9, 4/3, 3, 4]
print(f"{NAME}: n={n}, {S} samples")
for w in sorted(want):
    j = np.argmin(np.abs(u - w))
    if abs(u[j] - w) < 1e-5: print(f"  d^2={u[j]:.6f} pairs={cnt[j]:6d} P(same)={tot[j]/cnt[j]:.4f}")
order = np.argsort(-tot / np.maximum(cnt, 1))
print("  MOST often monochromatic (>=200 pairs):")
c = 0
for j in order:
    if cnt[j] >= 200 and u[j] > 1e-9:
        print(f"    d^2={u[j]:.6f} pairs={cnt[j]:6d} P(same)={tot[j]/cnt[j]:.4f}"); c += 1
        if c == 8: break
