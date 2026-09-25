"""Stack translates.  The first one moved the needle; stacking should multiply.

A single translate union of G costs 723211 conflicts to sweep with its dearest
pair at 3689, against G's own 6410 and 29 -- a hundredfold in both.  For scale,
Y's UNFORCED pairs at four colours cost 157 to 35365 conflicts and its forced
one over two million.  So one translate lifts G from nowhere near the band
where forcing lives into the bottom of it.

Stacking is the obvious next move and costs nothing arithmetically: G, G + t1,
G + t2, ... all at once, every pair of copies contributing its own cross
edges.  Report the conflict cost at each depth, because that is the quantity
that is actually moving.
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
          "39c59179b415/scratchpad/gtrans.pkl", "rb") as fh:
    TR = pickle.load(fh)
print(f"G: {len(P)} points; stacking translates, sampled counts "
      f"{[c for c, _ in TR[:6]]}  [{time.time()-t0:.0f}s]", flush=True)

for depth in (1, 2, 3, 4, 5):
    shifts = [((Fr(0),) * DIM, (Fr(0),) * DIM)] + [t for _, t in TR[:depth]]
    allp, allz, idx, copies = [], [], {}, []
    for sh in shifts:
        shz = (fl(sh[0]), fl(sh[1]))
        this = []
        for p, z in zip(P, zf):
            q = (madd(p[0], sh[0]), madd(p[1], sh[1]))
            if q in idx:
                this.append(idx[q])
            else:
                idx[q] = len(allp)
                this.append(len(allp))
                allp.append(q)
                allz.append((z[0] + shz[0], z[1] + shz[1]))
        copies.append(this)
    n = len(allp)
    E = set()
    for c in copies:
        for a0, b0 in GE:
            x, y = c[a0], c[b0]
            E.add((min(x, y), max(x, y)))
    inside = len(E)
    cell = {}
    for i in range(n):
        a, b = allz[i]
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    for i in range(n):
        a, b = allz[i]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i or (i, j) in E:
                        continue
                    if abs((a - allz[j][0]) ** 2
                           + (b - allz[j][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = (msub(allp[i][0], allp[j][0]),
                         msub(allp[i][1], allp[j][1]))
                    if norm2(d) == ONE:
                        E.add((i, j))
    E = sorted(E)
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** depth {depth}: {n} points, {len(E)} edges, "
              f"{len(E)-inside} cross -- NOT 5-COLOURABLE ***  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/tsix.pkl", "wb") as fh:
            pickle.dump((depth, allp, E), fh)
        break
    print(f"  depth {depth}: {n} points, {len(E)} edges "
          f"({len(E)-inside} cross), 5-colourable  [{time.time()-t0:.0f}s]",
          flush=True)
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
        sv.conf_budget(60000)
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
    print(f"    {len(cand)} pairs, {len(hits)} FORCED, {len(hard)} rerun, "
          f"{tot} conflicts, dearest {dear}  (G alone: 6410 and 29; Y's "
          f"unforced band: 157 to 35365)  [{time.time()-t0:.0f}s]",
          flush=True)
    if hits:
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/tforced.pkl", "wb") as fh:
            pickle.dump((depth, allp, E, hits), fh)
        print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
              flush=True)
        break
