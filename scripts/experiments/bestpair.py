"""Rank by what the pair already does, not by what the ring looks like.

The graded metric changes the picture.  Sa at four colours has a non-edge pair
agreeing in 84 per cent of sampled colourings against a chance of 25 -- 3.36
times over -- and rho pushes it to 100 per cent with six cross edges.  G at
five colours has a pair agreeing in 66 per cent against a chance of 20: 3.30
times over.  The same relative position, before the rotation.

Every ring scan here ranked rings by population, or by which radicals they
pay, and never by whether the pair they aim at already agrees.  That is the
one ordering the mechanism actually cares about.

So: find G's best-agreeing pairs at five colours; for each, take its midpoint
and the squared radius of the ring it sits on; and if that ring is rational
and closable, build G u rho(G) about the midpoint and see whether the pair's
agreement moves the way Sa's did.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SAMPLES, k = 32, 5
SC = ("/tmp/hn/")
t0 = time.time()
P = build_G(F, as_graph=False)
g = build_graph(P)
n = g.n
E = set((min(a, b), max(a, b)) for a, b in g.edges())
GE = sorted(E)
ONE = F.rational(1)


def colourings(pts, edges, kk, samples=SAMPLES, seed=2718):
    m = len(pts)
    cls = [[1 + v * kk + c for c in range(kk)] for v in range(m)]
    for a, b in edges:
        for c in range(kk):
            cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None, None
    rng, cols = random.Random(seed), []
    for s in range(samples):
        sv.set_phases([-(1 + w) if rng.random() < .05 else (1 + w)
                       for w in range(m * kk)])
        sv.solve()
        mo = sv.get_model()
        cols.append([next(c for c in range(kk) if mo[w * kk + c] > 0)
                     for w in range(m)])
    return np.array(cols, dtype=np.int8), sv


C, sv = colourings(P, GE, k)
sv.delete()
top = []
for i in range(n - 1):
    agree = (C[:, i + 1:] == C[:, i:i + 1]).sum(axis=0)
    for off in np.nonzero(agree >= 0.55 * SAMPLES)[0]:
        j = int(off) + i + 1
        if (i, j) not in E:
            top.append((int(agree[off]), i, j))
top.sort(reverse=True)
print(f"G at five: {len(top)} non-edge pairs agreeing in at least "
      f"{int(0.55*SAMPLES)}/{SAMPLES}; best {top[0][0] if top else 0}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

# For each, the ring it sits on about its own midpoint.
cand = []
for a, i, j in top[:60]:
    mid = Point((P[i].x + P[j].x) / F.rational(2),
                (P[i].y + P[j].y) / F.rational(2))
    d2 = (P[i].x - mid.x) ** 2 + (P[i].y - mid.y) ** 2
    v = float(d2)
    D = Fr(round(v * 55440), 55440)
    if abs(float(D) - v) > 1e-9 or D < Fr(1, 4):
        continue
    if closable_distance(D):
        cand.append((a, i, j, mid, D))
print(f"{len(cand)} of the top 60 sit on a rational closable ring about "
      f"their own midpoint: {[(a, str(D)) for a, _, _, _, D in cand[:8]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

for a, i, j, mid, D in cand[:8]:
    rho = rotation_joining(D, F).about(mid)
    U, idx = list(P), {p: q for q, p in enumerate(P)}
    for p in P:
        q = rho(p)
        if q not in idx:
            idx[q] = len(U)
            U.append(q)
    zf = [(float(p.x), float(p.y)) for p in U]
    ed = set(GE)
    img = [idx[rho(p)] for p in P]
    for x, y in GE:
        u, w = img[x], img[y]
        ed.add((min(u, w), max(u, w)))
    cell = {}
    for q, (x, y) in enumerate(zf):
        cell.setdefault((int(x // 1), int(y // 1)), []).append(q)
    cross = 0
    for q in range(n, len(U)):
        x, y = zf[q]
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for w in cell.get((int(x // 1) + dx, int(y // 1) + dy), ()):
                    if w >= q or abs((x - zf[w][0]) ** 2
                                     + (y - zf[w][1]) ** 2 - 1) > 1e-7:
                        continue
                    if ((U[w].x - U[q].x) ** 2 + (U[w].y - U[q].y) ** 2
                            == ONE) and (w, q) not in ed:
                        ed.add((w, q))
                        cross += 1
    CU, sv2 = colourings(U, sorted(ed), k)
    if CU is None:
        print(f"  *** pair {a}/{SAMPLES}, D={D}: union NOT 5-COLOURABLE ***",
              flush=True)
        break
    now = int((CU[:, i] == CU[:, j]).sum())
    forced = not sv2.solve(assumptions=[1 + i * k, -(1 + j * k)])
    sv2.delete()
    star = "  *** FORCED ***" if forced else ""
    print(f"  pair {i},{j}: {a}/{SAMPLES} in G -> {now}/{SAMPLES} in the "
          f"union (D={D}, {len(U)} points, {cross} cross){star}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if forced:
        with open(SC + "bestpair_forced.pkl", "wb") as fh:
            pickle.dump((str(D), i, j,
                         [(str(p.x), str(p.y)) for p in U], sorted(ed)), fh)
        break
