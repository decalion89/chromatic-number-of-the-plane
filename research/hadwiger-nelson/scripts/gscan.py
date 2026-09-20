"""G u u.G for the rotations that actually bite, at five colours.

The four-colour step was a union whose rotation bites.  This is the same move
one floor up, with the rotations found by solving rather than by guessing
which distances G realises -- the earlier attempt used closing rotations and
got +5 cross edges on a ring of eleven, which is not the experiment.

Two phases.  First the cheap probe on every union: how much does it really
bite, and is it still 5-colourable?  A union that is not is the answer to the
whole problem.  Then the pair scan on the ones that bite hardest, with the
protocol the budget lesson forced -- collect the undecided, rerun them with no
budget at all, and only then read the result.
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from msqrt import madd, msub, mmul
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.homcol import closable_distance
from pysat.solvers import Solver

GENS = (3, 5, 7, 11)
DIM = 16
t0 = time.time()
pts = build_G(F, as_graph=False)
P = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c)) for p in pts]
zf = [(float(p.x), float(p.y)) for p in pts]
ONE = (Fr(1),) + (Fr(0),) * (DIM - 1)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/gsolve.pkl", "rb") as fh:
    ROT = pickle.load(fh)
print(f"G: {len(P)} points, {len(ROT)} biting rotations  "
      f"[{time.time()-t0:.0f}s]", flush=True)
RT = [[float(sum(float(c) * (r ** 0.5 if r > 1 else 1)
                 for c, r in zip(comp, RAD))) for comp in u]
      for u in ROT] if False else None
RAD = []
for m in range(DIM):
    pr = 1
    for b in range(4):
        if m >> b & 1:
            pr *= GENS[b]
    RAD.append(pr ** 0.5)


def fl(c):
    return sum(float(x) * r for x, r in zip(c, RAD))


def cmulf(u, z):
    return (fl(u[0]) * z[0] - fl(u[1]) * z[1],
            fl(u[1]) * z[0] + fl(u[0]) * z[1])


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


def cmul(z, w):
    return (msub(mmul(z[0], w[0], GENS), mmul(z[1], w[1], GENS)),
            madd(mmul(z[0], w[1], GENS), mmul(z[1], w[0], GENS)))


rows = []
for ri, u in enumerate(ROT):
    img = [cmul(u, p) for p in P]
    zi = [cmulf(u, z) for z in zf]
    allp, allz, seen = list(P), list(zf), set(P)
    shared = 0
    for q, z in zip(img, zi):
        if q in seen:
            shared += 1
        else:
            seen.add(q)
            allp.append(q)
            allz.append(z)
    cell = {}
    for i, (a, b) in enumerate(allz):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    E, cross = [], 0
    n0 = len(P)
    for i, (a, b) in enumerate(allz):
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i:
                        continue
                    if abs((a - allz[j][0]) ** 2
                           + (b - allz[j][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = (msub(allp[i][0], allp[j][0]),
                         msub(allp[i][1], allp[j][1]))
                    if norm2(d) == ONE:
                        E.append((i, j))
                        if (i < n0) != (j < n0):
                            cross += 1
    cls = [[1 + v * 5 + c for c in range(5)] for v in range(len(allp))]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        five = sv.solve()
    rows.append((cross, shared, len(allp), len(E), ri))
    if not five:
        print(f"  *** rotation {ri}: {len(allp)} points, {len(E)} edges, "
              f"{cross} cross -- NOT 5-COLOURABLE ***  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/gsix.pkl", "wb") as fh:
            pickle.dump((ri, allp, E), fh)
        break
    if ri % 10 == 0:
        best = max(rows)[0] if rows else 0
        print(f"  ... {ri}/{len(ROT)}, best bite so far {best} cross edges "
              f"(last: {cross} cross, {shared} shared, {len(allp)} pts)  "
              f"[{time.time()-t0:.0f}s]", flush=True)
rows.sort(reverse=True)
print(f"\ncross-edge counts, best first: {[r[0] for r in rows[:15]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/gscan.pkl", "wb") as fh:
    pickle.dump(rows, fh)
