"""The same architecture, one floor up.

de Grey's tower is: H, the hexagon with its centre; Sa, its closure under the
12-element dihedral group; Y = Sa u rho(Sa), which forces a pair at four
colours; G = Y u rho'(Y), which is 5-chromatic.  Each floor is built from the
floor below by closing under rotations and then taking one union.

So the floor above starts from G itself.  W is G closed under the same
dihedral group -- twelve images about the origin -- and the question is
whether that density is enough to force a pair at FIVE colours, which G alone
does not (21358 pairs tested, none forced, none even hard).

The pair test is the same: fixing c(a) = 0 is free, so each query is one pair
of assumptions against an incremental solver, and unsatisfiable is the prize.
A budget flags the hard ones instead of letting a single query stall the scan.
"""
import sys, time
from fractions import Fraction as Fr
from collections import Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()
G = build_G(F, as_graph=False)
rot60 = _rot60(F)
W, seen = [], set()
for p in G:
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                W.append(q)
            q = rot60(q)
print(f"W: {len(W)} points from G's {len(G)}  [{time.time()-t0:.0f}s]",
      flush=True)

zs = [(float(p.x), float(p.y)) for p in W]
cell = {}
for i, (a, b) in enumerate(zs):
    cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
cands = []
for i, (a, b) in enumerate(zs):
    cx, cy = int(a // 1), int(b // 1)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            for j in cell.get((cx + da, cy + db), ()):
                if j > i and abs((a - zs[j][0]) ** 2
                                 + (b - zs[j][1]) ** 2 - 1) < 1e-7:
                    cands.append((i, j))
print(f"  {len(cands)} candidate edges  [{time.time()-t0:.0f}s]", flush=True)
E = []
for i, j in cands:
    d = W[i] - W[j]
    if d.x * d.x + d.y * d.y == F.rational(1):
        E.append((i, j))
print(f"  {len(E)} edges, mean degree {2*len(E)/len(W):.1f} "
      f"(G's was 10.0)  [{time.time()-t0:.0f}s]", flush=True)

K = 5
cls = [[1 + v * K + c for c in range(K)] for v in range(len(W))]
for a, b in E:
    for c in range(K):
        cls.append([-(1 + a * K + c), -(1 + b * K + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
ok = sv.solve()
print(f"  5-colourable? {ok}  [{time.time()-t0:.0f}s]", flush=True)
if not ok:
    print("  *** W IS 6-CHROMATIC ***", flush=True)
    sys.exit(0)

deg = Counter()
for a, b in E:
    deg[a] += 1
    deg[b] += 1
hub = [v for v, _ in deg.most_common(400)]
print(f"  scanning pairs among the 400 busiest vertices (degree "
      f"{deg[hub[0]]} down to {deg[hub[-1]]})  [{time.time()-t0:.0f}s]",
      flush=True)
cand2 = []
for x in range(len(hub)):
    i = hub[x]
    ai, bi = zs[i]
    for y in range(x + 1, len(hub)):
        j = hub[y]
        v = (ai - zs[j][0]) ** 2 + (bi - zs[j][1]) ** 2
        if v > 36.0:
            continue
        D = Fr(round(v * 1584), 1584)
        if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
            continue
        cand2.append((i, j, D))
print(f"  {len(cand2)} pairs at a closable distance  "
      f"[{time.time()-t0:.0f}s]", flush=True)
hits, hard = [], []
for k, (i, j, D) in enumerate(cand2):
    sv.conf_budget(30000)
    r = sv.solve_limited(assumptions=[1 + i * K, -(1 + j * K)])
    if r is False:
        hits.append((i, j, D))
        print(f"  *** FORCED AT FIVE COLOURS: {i},{j} D = {D}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    elif r is None:
        hard.append((i, j, D))
    if k and k % 2000 == 0:
        print(f"  ... {k}/{len(cand2)}, {len(hard)} hard, {len(hits)} forced "
              f" [{time.time()-t0:.0f}s]", flush=True)
for i, j, D in hard:
    if not sv.solve(assumptions=[1 + i * K, -(1 + j * K)]):
        hits.append((i, j, D))
        print(f"  *** FORCED (full run): {i},{j} D = {D}", flush=True)
print(f"  {len(hits)} forced pairs at five colours in W ({len(hard)} needed "
      f"a full run)  [{time.time()-t0:.0f}s]", flush=True)
