"""The pair scan on unions that really bite, at five colours.

Rotating G about a vertex of degree 60 instead of about the origin changes
everything the solve returns: 796 and 1182 cross edges against the one or two
the origin gives, and against the six that make de Grey's own step work.  All
of them still 5-colourable -- which is what a candidate looks like, not what a
failure looks like.

So scan them.  A pair forced monochromatic in every 5-colouring, at a distance
the field closes, is a 6-chromatic unit-distance graph: spindle it and both
copies read the same colour on two adjacent points.

Protocol as the budget lesson forced it -- collect every undecided query and
rerun it with no budget at all before reading anything.
"""
import sys, time, pickle, glob
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
P0 = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c))
      for p in pts]
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


def cmul(z, w):
    return (msub(mmul(z[0], w[0], GENS), mmul(z[1], w[1], GENS)),
            madd(mmul(z[0], w[1], GENS), mmul(z[1], w[0], GENS)))


for path in sorted(glob.glob("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-"
                             "f432-5848-a506-39c59179b415/scratchpad/"
                             "gpivot_*.pkl")):
    with open(path, "rb") as fh:
        pv, ROT = pickle.load(fh)
    o = P0[pv]
    P = [(msub(p[0], o[0]), msub(p[1], o[1])) for p in P0]
    oz = (fl(o[0]), fl(o[1]))
    zt = [(z[0] - oz[0], z[1] - oz[1]) for z in zf]
    print(f"\npivot {pv}: {len(ROT)} rotations  [{time.time()-t0:.0f}s]",
          flush=True)
    for ri, u in enumerate(ROT):
        ux, uy = fl(u[0]), fl(u[1])
        img = [cmul(u, p) for p in P]
        zi = [(ux * z[0] - uy * z[1], uy * z[0] + ux * z[1]) for z in zt]
        allp, allz, seen = list(P), list(zt), set(P)
        idx = {p: i for i, p in enumerate(P)}
        second = []
        for q, z in zip(img, zi):
            if q in seen:
                second.append(idx[q])
            else:
                seen.add(q)
                idx[q] = len(allp)
                second.append(len(allp))
                allp.append(q)
                allz.append(z)
        n0 = len(P)
        cell = {}
        for i in range(len(allz)):
            a, b = allz[i]
            cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
        E = set()
        for a0, b0 in GE:
            E.add((min(a0, b0), max(a0, b0)))
            x, y = second[a0], second[b0]
            E.add((min(x, y), max(x, y)))
        cross = 0
        for j in range(n0, len(allz)):
            a, b = allz[j]
            cx, cy = int(a // 1), int(b // 1)
            for da in (-1, 0, 1):
                for db in (-1, 0, 1):
                    for i in cell.get((cx + da, cy + db), ()):
                        if i >= n0 or abs((a - allz[i][0]) ** 2
                                          + (b - allz[i][1]) ** 2 - 1) > 1e-7:
                            continue
                        d = (msub(allp[i][0], allp[j][0]),
                             msub(allp[i][1], allp[j][1]))
                        if norm2(d) == ONE:
                            E.add((i, j))
                            cross += 1
        E = sorted(E)
        n = len(allp)
        cls = [[1 + v * 5 + c for c in range(5)] for v in range(n)]
        for a, b in E:
            for c in range(5):
                cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        if not sv.solve():
            print(f"  *** rotation {ri}: {n} points, {len(E)} edges NOT "
                  f"5-COLOURABLE ***  [{time.time()-t0:.0f}s]", flush=True)
            with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                      "a506-39c59179b415/scratchpad/gsix.pkl", "wb") as fh:
                pickle.dump((pv, ri, allp, E), fh)
            sys.exit(0)
        cand = []
        for i in range(n):
            ai, bi = allz[i]
            for j in range(i + 1, n):
                v = (ai - allz[j][0]) ** 2 + (bi - allz[j][1]) ** 2
                if v > 36.0:
                    continue
                D = Fr(round(v * 1584), 1584)
                if abs(float(D) - v) > 1e-7 or D == 1 \
                        or not closable_distance(D):
                    continue
                cand.append((i, j))
        hits, hard = [], []
        for i, j in cand:
            sv.conf_budget(30000)
            r = sv.solve_limited(assumptions=[1 + i * 5, -(1 + j * 5)])
            if r is False:
                hits.append((i, j))
            elif r is None:
                hard.append((i, j))
        for i, j in hard:
            if not sv.solve(assumptions=[1 + i * 5, -(1 + j * 5)]):
                hits.append((i, j))
        sv.delete()
        print(f"  rotation {ri}: {n} pts, {len(E)} edges, {cross} cross, "
              f"{len(cand)} pairs, {len(hits)} FORCED, {len(hard)} rerun  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if hits:
            with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                      "a506-39c59179b415/scratchpad/gforced.pkl", "wb") as fh:
                pickle.dump((pv, ri, allp, E, hits), fh)
            print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
                  flush=True)
            sys.exit(0)
