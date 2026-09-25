"""Stack the census leaders, and measure the decay profile at every depth.

Conflict counts grow with the graph, so they compare nothing across sizes.
The survivor curve does.  Sample proper 5-colourings one at a time and watch
how many candidate pairs still agree: after k samples the count says how
tightly the union constrains its own colourings, in units that are pairs, not
solver internals.  A union near forcing decays slowly and never reaches zero,
because a forced pair survives every sample by definition.  A union far from
forcing collapses in a dozen.

Y is the calibration and it cannot be argued with: Y has a forced pair, so
Y's curve has a floor at 1.  Anything whose curve reaches 0 has no forced
pair at all, witnessed, and that is a complete negative rather than an
exhausted budget.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import random
from msqrt import madd, msub, mmul
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

GENS = (3, 5, 7, 11)
DIM = 16
DEPTH = 8
MARKS = (1, 2, 3, 5, 10, 20, 40)
SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
random.seed(11235)
t0 = time.time()
pts = build_G(F, as_graph=False)
P = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c)) for p in pts]
zf = [(float(p.x), float(p.y)) for p in pts]
GE = list(build_graph(pts).edges())
ONE = (Fr(1),) + (Fr(0),) * (DIM - 1)
RAD = []
for m in range(DIM):
    pr = 1
    for b in range(4):
        if m >> b & 1:
            pr *= GENS[b]
    RAD.append(pr ** 0.5)


def fl(c):
    return sum(float(x) * r for x, r in zip(c, RAD))


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


def profile(sv, n, cand, marks=MARKS, cap=60, stall=20):
    """Sample colourings; return the survivor curve and the final survivors."""
    surv, curve, since = cand, [], 0
    for s in range(cap):
        for _t in range(6):
            vs = random.sample(range(n), 8)
            if sv.solve(assumptions=[1 + v * 5 + random.randrange(5)
                                     for v in vs]):
                break
        else:
            sv.solve()
        m = sv.get_model()
        col = [0] * n
        for w in range(n):
            for c in range(5):
                if m[w * 5 + c] > 0:
                    col[w] = c
                    break
        before = len(surv)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        since = 0 if len(surv) < before else since + 1
        if s + 1 in marks:
            curve.append((s + 1, len(surv)))
        if not surv or since >= stall:
            break
    curve.append((s + 1, len(surv)))
    return curve, surv


with open(SC + "tcensus.pkl", "rb") as fh:
    CEN = pickle.load(fh)
print(f"census leaders: {[c for c, _ in CEN[:DEPTH]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

allp, allz = list(P), list(zf)
idx = {p: i for i, p in enumerate(P)}
E = set((min(a, b), max(a, b)) for a, b in GE)
for depth in range(1, DEPTH + 1):
    t = CEN[depth - 1][1]
    tz = (fl(t[0]), fl(t[1]))
    base = len(allp)
    new = []
    for p, z in zip(P, zf):
        q = (madd(p[0], t[0]), madd(p[1], t[1]))
        if q in idx:
            new.append(idx[q])
        else:
            idx[q] = len(allp)
            new.append(len(allp))
            allp.append(q)
            allz.append((z[0] + tz[0], z[1] + tz[1]))
    for a0, b0 in GE:
        x, y = new[a0], new[b0]
        E.add((min(x, y), max(x, y)))
    cell = {}
    for i in range(len(allp)):
        a, b = allz[i]
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    cross = 0
    for j in range(base, len(allp)):
        a, b = allz[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if i >= j or abs((a - allz[i][0]) ** 2
                                     + (b - allz[i][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = (msub(allp[i][0], allp[j][0]),
                         msub(allp[i][1], allp[j][1]))
                    if norm2(d) == ONE and (i, j) not in E:
                        E.add((i, j))
                        cross += 1
    n = len(allp)
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in sorted(E):
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    print(f"  depth {depth}: {n} points, {len(E)} edges (+{cross} cross), "
          f"5-colourable {ok}  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("  *** NOT 5-COLOURABLE -- SIX COLOURS ***", flush=True)
        with open(SC + "tsfast_six.pkl", "wb") as fh:
            pickle.dump((depth, allp, sorted(E)), fh)
        sys.exit(0)
    cand = []
    for i in range(n):
        ai, bi = allz[i]
        for j in range(i + 1, n):
            v = (ai - allz[j][0]) ** 2 + (bi - allz[j][1]) ** 2
            if v > 36.0:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
                continue
            cand.append((i, j))
    curve, surv = profile(sv, n, cand)
    hits = [(i, j) for i, j in surv
            if not sv.solve(assumptions=[1 + i * 5, -(1 + j * 5)])]
    sv.delete()
    print(f"    {len(cand)} pairs, decay {curve}, {len(hits)} FORCED"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if hits:
        with open(SC + "tsfast_forced.pkl", "wb") as fh:
            pickle.dump((depth, allp, sorted(E), hits), fh)
        print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
              flush=True)
        break
