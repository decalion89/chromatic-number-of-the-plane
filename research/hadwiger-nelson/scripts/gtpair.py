"""The translate unions, measured exactly and scanned at five colours.

The histogram of t = p - v - q over G puts the best non-trivial translate at
451 cross edges from a sample of 150 of the 1581 points, so the true count
should be near five thousand -- against the 1578 of the best rotation stack
and the six that make de Grey's own step work.

Translates cost nothing arithmetically: no square has to lie in the field, no
rotation has to exist.  They were simply never tried.

Build each union exactly, count what really crosses, and scan every pair at a
closable distance, rerunning anything the budget leaves undecided.
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
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
          "39c59179b415/scratchpad/gtrans.pkl", "rb") as fh:
    TR = pickle.load(fh)
print(f"G: {len(P)} points; {len(TR)} translates, sampled counts "
      f"{[c for c, _ in TR[:8]]}  [{time.time()-t0:.0f}s]", flush=True)

Pidx = {p: i for i, p in enumerate(P)}
for ti, (sc, t) in enumerate(TR):
    tz = (fl(t[0]), fl(t[1]))
    img = [(madd(p[0], t[0]), madd(p[1], t[1])) for p in P]
    allp, allz = list(P), list(zf)
    idx = dict(Pidx)
    second = []
    for q, z in zip(img, zf):
        if q in idx:
            second.append(idx[q])
        else:
            idx[q] = len(allp)
            second.append(len(allp))
            allp.append(q)
            allz.append((z[0] + tz[0], z[1] + tz[1]))
    n0, n = len(P), len(allp)
    shared = 2 * n0 - n
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
        print(f"  *** translate {ti}: {n} points, {len(E)} edges, {cross} "
              f"cross -- NOT 5-COLOURABLE ***  [{time.time()-t0:.0f}s]",
              flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/gsix_t.pkl", "wb") as fh:
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

    base, dear, hits, hard = conf(), 0, [], []
    for i, j in cand:
        sv.conf_budget(30000)
        r = sv.solve_limited(assumptions=[1 + i * 5, -(1 + j * 5)])
        if r is False:
            hits.append((i, j))
        elif r is None:
            hard.append((i, j))
        c2 = conf()
        dear = max(dear, c2 - base)
        base = c2
    for i, j in hard:
        if not sv.solve(assumptions=[1 + i * 5, -(1 + j * 5)]):
            hits.append((i, j))
    tot = conf()
    sv.delete()
    print(f"  translate {ti}: {n} pts ({shared} shared), {len(E)} edges, "
          f"{cross} cross, {len(cand)} pairs, {len(hits)} FORCED, "
          f"{len(hard)} rerun, {tot} conflicts, dearest {dear}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if hits:
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/gforced_t.pkl", "wb") as fh:
            pickle.dump((ti, allp, E, hits), fh)
        print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
              flush=True)
        break
