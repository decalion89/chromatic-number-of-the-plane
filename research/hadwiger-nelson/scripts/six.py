"""The same lemma, one floor up.

  If H carries a pair p, q at distance d that is monochromatic in every proper
  k-colouring, and the field admits a rotation rho about p closing d, then
  H u rho_p(H) needs k+1 colours.

At k = 4 that is de Grey.  At k = 5 it is chi(R^2) >= 6, and it asks a
question about a graph I already have: does his G force a pair at five
colours?

The distances that can be closed over his field Q(sqrt-3, sqrt-7, sqrt-11,
sqrt-15) are those with 1 - 4D in one of its sixteen rational square classes,
the squarefree products of subsets of {-3, -7, -11, -15}.  There are many
more of them than K has, which is the one advantage that field holds.

Most queries are satisfiable and cheap; a pair that is really forced makes
the solver prove unsatisfiability, so hard ones are flagged by a conflict
budget rather than run to the end, and revisited without a budget.
"""
import sys, time
from fractions import Fraction as Fr
from itertools import combinations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD
from pysat.solvers import Solver

t0 = time.time()
F = DEGREY_FIELD
G = build_G(F, as_graph=False)
n = len(G)
zs = [(float(p.x), float(p.y)) for p in G]
print(f"G: {n} vertices  [{time.time()-t0:.0f}s]", flush=True)


def squarefree(m):
    if m == 0:
        return 0
    s, d = (-1 if m < 0 else 1), abs(m)
    f = 2
    while f * f <= d:
        e = 0
        while d % f == 0:
            d //= f
            e += 1
        if e % 2:
            s *= f
        f += 1
    return s * d


CLASSES = set()
for r in range(5):
    for sub in combinations((-3, -7, -11, -15), r):
        pr = 1
        for v in sub:
            pr *= v
        CLASSES.add(squarefree(pr))
print(f"  square classes of de Grey's field: {sorted(CLASSES)}", flush=True)


def closable(D):
    r = Fr(1) - 4 * Fr(D)
    if r == 0:
        return False
    return squarefree(r.numerator * r.denominator) in CLASSES


cell = {}
for i, (a, b) in enumerate(zs):
    cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
E = []
for i, (a, b) in enumerate(zs):
    cx, cy = int(a // 1), int(b // 1)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            for j in cell.get((cx + da, cy + db), ()):
                if j <= i:
                    continue
                if abs((a - zs[j][0]) ** 2 + (b - zs[j][1]) ** 2 - 1) > 1e-7:
                    continue
                d = G[i] - G[j]
                if d.x * d.x + d.y * d.y == F.rational(1):
                    E.append((i, j))
print(f"  {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

K = 5
cls = [[1 + v * K + c for c in range(K)] for v in range(n)]
for a, b in E:
    for c in range(K):
        cls.append([-(1 + a * K + c), -(1 + b * K + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
print(f"  5-colourable? {sv.solve()}  [{time.time()-t0:.0f}s]", flush=True)

cand = []
for i in range(n):
    ai, bi = zs[i]
    for j in range(i + 1, n):
        v = (ai - zs[j][0]) ** 2 + (bi - zs[j][1]) ** 2
        if v > 20.0:
            continue
        D = Fr(round(v * 1584), 1584)
        if abs(float(D) - v) > 1e-7 or D == 1 or not closable(D):
            continue
        cand.append((i, j, D))
print(f"  {len(cand)} pairs at a closable distance  "
      f"[{time.time()-t0:.0f}s]", flush=True)

hard, hits = [], []
for k, (i, j, D) in enumerate(cand):
    sv.conf_budget(30000)
    r = sv.solve_limited(assumptions=[1 + i * K, -(1 + j * K)])
    if r is False:
        hits.append((i, j, D))
        print(f"  *** FORCED AT FIVE COLOURS: {i},{j} D = {D}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    elif r is None:
        hard.append((i, j, D))
    if k and k % 2000 == 0:
        print(f"  ... {k}/{len(cand)}, {len(hard)} over budget, "
              f"{len(hits)} forced  [{time.time()-t0:.0f}s]", flush=True)
print(f"  pass 1: {len(hits)} forced, {len(hard)} over budget  "
      f"[{time.time()-t0:.0f}s]", flush=True)
for i, j, D in hard:
    if not sv.solve(assumptions=[1 + i * K, -(1 + j * K)]):
        hits.append((i, j, D))
        print(f"  *** FORCED AT FIVE COLOURS (full run): {i},{j} D = {D}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"  total {len(hits)} forced pairs at five colours  "
      f"[{time.time()-t0:.0f}s]", flush=True)
