"""Scan the translates in census order -- the heaviest crossing first.

The sampled histogram ranks them badly: its top four give 1442, 1440, 1442 and
1444 cross edges exactly, while one it ranked lower gives 1528, and that one
costs 17555 conflicts on its dearest pair against the others' 3689.  Five times
the cost for six percent more crossing.  So the census, which counts exactly,
is the ordering that matters.

Protocol as the budget lesson forced it, but used properly this time.  A
forced pair is UNSAT and will blow ANY budget, so a high budget on the sweep
costs nothing in correctness: everything it leaves undecided is reruns without
one, and only then does a verdict mean anything.  With the dearest observed at
17555, a budget of 200000 leaves almost nothing to rerun and turns a six-hour
sweep into minutes.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from msqrt import madd, msub, mmul
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

GENS = (3, 5, 7, 11)
DIM = 16
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


with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/tcensus.pkl", "rb") as fh:
    CEN = pickle.load(fh)
print(f"census: {len(CEN)} translates, cross counts "
      f"{[c for c, _ in CEN[:12]]}  [{time.time()-t0:.0f}s]", flush=True)

Pidx = {p: i for i, p in enumerate(P)}
for ti, (ccount, t) in enumerate(CEN[:24]):
    tz = (fl(t[0]), fl(t[1]))
    allp, allz = list(P), list(zf)
    idx = dict(Pidx)
    second = []
    for p, z in zip(P, zf):
        q = (madd(p[0], t[0]), madd(p[1], t[1]))
        if q in idx:
            second.append(idx[q])
        else:
            idx[q] = len(allp)
            second.append(len(allp))
            allp.append(q)
            allz.append((z[0] + tz[0], z[1] + tz[1]))
    n0, n = len(P), len(allp)
    cell = {}
    for i in range(n0):
        a, b = allz[i]
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    E = set()
    for a0, b0 in GE:
        E.add((min(a0, b0), max(a0, b0)))
        x, y = second[a0], second[b0]
        E.add((min(x, y), max(x, y)))
    cross = 0
    for j in range(n0, n):
        a, b = allz[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if abs((a - allz[i][0]) ** 2
                           + (b - allz[i][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = (msub(allp[i][0], allp[j][0]),
                         msub(allp[i][1], allp[j][1]))
                    if norm2(d) == ONE and (i, j) not in E:
                        E.add((i, j))
                        cross += 1
    E = sorted(E)
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** translate {ti}: {n} points, {cross} cross -- NOT "
              f"5-COLOURABLE ***  [{time.time()-t0:.0f}s]", flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/tbest_six.pkl", "wb") as fh:
            pickle.dump((ti, allp, E), fh)
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

    def conf():
        return sv.accum_stats().get("conflicts", 0)

    base, dear, hard = conf(), 0, []
    for i, j in cand:
        sv.conf_budget(200000)
        r = sv.solve_limited(assumptions=[1 + i * 5, -(1 + j * 5)])
        if r is None:
            hard.append((i, j))
        c2 = conf()
        dear = max(dear, c2 - base)
        base = c2
    hits = []
    for i, j in hard:
        if not sv.solve(assumptions=[1 + i * 5, -(1 + j * 5)]):
            hits.append((i, j))
    tot = conf()
    sv.delete()
    print(f"  translate {ti}: {n} pts, {len(E)} edges, {cross} cross, "
          f"{len(cand)} pairs, {len(hits)} FORCED, {len(hard)} rerun, "
          f"{tot} conflicts, dearest {dear}  [{time.time()-t0:.0f}s]",
          flush=True)
    if hits:
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/tbest_forced.pkl", "wb") as fh:
            pickle.dump((ti, allp, E, hits), fh)
        print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
              flush=True)
        break
