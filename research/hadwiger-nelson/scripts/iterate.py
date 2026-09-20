"""Iterate the recursion: make Z the new core and ask its rings.

Z = G* u rho_16(G*) is 20809 points, contains G so it is 5-chromatic, and is
invariant under the 60-degree turn about (-2,0) -- rho_16 shares that pivot
and so commutes with it.  That is the same shape G* had, one step further
along, and the criterion applies to it unchanged.

Building the spindle is no longer the cheap way to ask.  If the union forces
its ring's antipodal pair, the spindle is uncolourable; if it does not, the
spindle colours -- so the sampling answers the same question without ever
building a graph of eighty thousand points.  A pair separated by a sampled
colouring is not forced, witnessed, and only a pair that survives every sample
is worth the spindle.
"""
import sys, time, pickle, random
from collections import Counter
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.field import Field, embed
from hn.geometry import DEGREY_FIELD as K, Point, _rot60, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
K17 = Field((3, 5, 7, 11, 17))
RAD17 = (3, 5, 7, 11, 17)
t0 = time.time()
PIV = Point(K17.rational(-2), K17.zero())

rot = _rot60(K).about(Point(K.rational(-2), K.zero()))
G = build_G(K, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
Gs = [Point(embed(p.x, K17), embed(p.y, K17)) for p in Gs]
rho16 = rotation_joining(16, K17).about(PIV)
Z, zs = list(Gs), set(Gs)
for p in Gs:
    q = rho16(p)
    if q not in zs:
        zs.add(q)
        Z.append(q)
print(f"Z: {len(Z)} points  [{time.time()-t0:.0f}s]", flush=True)

rings = Counter()
for p in Z:
    dx, dy = p.x - PIV.x, p.y - PIV.y
    rings[float(dx * dx + dy * dy)] += 1
cand = []
for r2f, cnt in rings.items():
    D = Fr(round(r2f * 720720), 720720)
    if abs(float(D) - r2f) > 1e-9 or D < Fr(1, 4) or D == 1:
        continue
    if closable_distance(D, RAD17) and closable_distance(4 * D, RAD17):
        cand.append((cnt, D))
cand.sort(reverse=True)
print(f"{len(rings)} rings of Z, {len(cand)} doubly usable over K(sqrt17): "
      f"{[(c, str(d)) for c, d in cand]}  [{time.time()-t0:.0f}s]",
      flush=True)

for cnt, D in cand[:6]:
    rho = rotation_joining(D, K17).about(PIV)
    U, us = list(Z), set(Z)
    for p in Z:
        q = rho(p)
        if q not in us:
            us.add(q)
            U.append(q)
    idx = {p: i for i, p in enumerate(U)}
    anti = []
    for i, p in enumerate(U):
        if (p.x - PIV.x) ** 2 + (p.y - PIV.y) ** 2 != K17.rational(D):
            continue
        j = idx.get(Point(PIV.x + PIV.x - p.x, PIV.y + PIV.y - p.y))
        if j is not None and i < j:
            anti.append((i, j))
    n = len(U)
    E = sorted((min(a, b), max(a, b)) for a, b in build_graph(U).edges())
    k = 5
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  D={D}: {n} points, {len(E)} edges -- NOT 5-COLOURABLE",
              flush=True)
        with open(SC + "iterate_six.pkl", "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in U], E), fh)
        sv.delete()
        break
    rng, surv, hist = random.Random(4649), list(anti), []
    for s in range(14):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + v)
                       for v in range(n * k)])
        sv.solve()
        m = sv.get_model()
        col = {}
        for i, j in surv:
            for w in (i, j):
                if w not in col:
                    col[w] = next(c for c in range(k) if m[w * k + c] > 0)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        hist.append(len(surv))
        if not surv:
            break
    sv.delete()
    print(f"  D={D} ({cnt} on the ring): {n} points, {len(E)} edges, "
          f"{len(anti)} antipodal, curve {hist}  [{time.time()-t0:.0f}s]",
          flush=True)
    if surv:
        print(f"  *** {len(surv)} antipodal pairs survive -- worth the "
              f"spindle ***", flush=True)
        with open(SC + f"iterate_surv_{D.numerator}_{D.denominator}.pkl",
                  "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in U], E, surv),
                        fh)
