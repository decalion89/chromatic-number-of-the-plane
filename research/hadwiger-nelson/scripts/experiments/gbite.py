"""Measure the pivot rotations exactly: points, shared, cross edges, colours.

About the origin a rotation of de Grey's field bites G with one or two cross
edges.  About a vertex of degree 60 the same solve turns up rotations serving
199 of 40 sampled q-values -- and 40 is one fortieth of the point set, so the
true count should be in the thousands.  His own step bites with six.

So build the union and count it exactly: how many points, how many shared, how
many edges really cross, and whether five colours still suffice.  A union that
needs six is the whole problem.
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


import glob
for path in sorted(glob.glob("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-"
                             "f432-5848-a506-39c59179b415/scratchpad/"
                             "gpivot_*.pkl")):
    with open(path, "rb") as fh:
        pv, ROT = pickle.load(fh)
    o = P0[pv]
    P = [(msub(p[0], o[0]), msub(p[1], o[1])) for p in P0]
    oz = (fl(o[0]), fl(o[1]))
    zt = [(z[0] - oz[0], z[1] - oz[1]) for z in zf]
    print(f"\npivot {pv}: {len(ROT)} candidate rotations  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for ri, u in enumerate(ROT[:6]):
        ux, uy = fl(u[0]), fl(u[1])
        img = [cmul(u, p) for p in P]
        zi = [(ux * z[0] - uy * z[1], uy * z[0] + ux * z[1]) for z in zt]
        allp, allz, seen = list(P), list(zt), set(P)
        shared = 0
        for q, z in zip(img, zi):
            if q in seen:
                shared += 1
            else:
                seen.add(q)
                allp.append(q)
                allz.append(z)
        n0 = len(P)
        cell = {}
        for i in range(n0):
            a, b = allz[i]
            cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
        crossE = []
        for j in range(n0, len(allz)):
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
                        if norm2(d) == ONE:
                            crossE.append((i, j))
        E = list(GE) + [(i + n0, j + n0) for i, j in GE if i + n0 < len(allp)
                        and j + n0 < len(allp)] + crossE
        cls = [[1 + v * 5 + c for c in range(5)] for v in range(len(allp))]
        for a, b in E:
            for c in range(5):
                cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            five = sv.solve()
        print(f"  rotation {ri}: {len(allp)} points ({shared} shared), "
              f"{len(crossE)} cross edges, 5-colourable {five}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if not five:
            print("  *** SIX COLOURS ***", flush=True)
            with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                      "a506-39c59179b415/scratchpad/gsix.pkl", "wb") as fh:
                pickle.dump((pv, ri, allp, E), fh)
            sys.exit(0)
