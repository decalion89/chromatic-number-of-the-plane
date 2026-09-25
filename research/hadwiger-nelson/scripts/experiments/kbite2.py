"""The 81 rotations that really bite, scanned for a forced pair.

Sampling found none -- all 3030 units of K leave Sa and its image sharing only
the origin, with zero cross edges once the float filter is checked exactly.
Solving finds 81.  For each pair p, q of Sa the rotation carrying q to
distance 1 from p is the root of an explicit quadratic, and it lies in K
exactly when (R^2 - A.P)/(-3) is a square in F = Q(m, sqrt33); 26136 pairs
give 81 rotations up to the zeta_6 symmetry Sa already has.

de Grey's own union has one shared point and six cross edges, and that is the
whole of what turns a graph with no forced pair into one that has one.  So:
build each union, count what really crosses, and scan every pair at a
K-closable distance.  A hit is a 5-chromatic unit-distance graph over a field
that blocks at the gate.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.degrey import S_POINTS
from hn.homcol import closable_over
from pysat.solvers import Solver

t0 = time.time()
KCL = (1, -3, -11, 33)
SQ3 = ksub(kmul(krat(2), Z6), K1)


def to_k(xs, ys):
    p, q = Fr(xs.get(1, 0)), Fr(xs.get(33, 0))
    r, t = Fr(ys.get(3, 0)), Fr(ys.get(11, 0))
    return kadd(kadd(krat(p), kmul(krat(r), SQ3)),
                kmul(ksub(krat(t), kmul(krat(q), SQ3)), S))


Sa, seenp = [], set()
for z in [to_k(x, y) for x, y in S_POINTS]:
    for base in (z, kconj(z)):
        w = base
        for _ in range(6):
            if w not in seenp:
                seenp.add(w)
                Sa.append(w)
            w = kmul(Z6, w)
saset = set(Sa)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/ksolve.pkl", "rb") as fh:
    raw = pickle.load(fh)
ROT = [tuple(tuple(v[3 * i:3 * i + 3]) for i in range(4)) for v in raw]
print(f"Sa {len(Sa)} points, {len(ROT)} solved rotations  "
      f"[{time.time()-t0:.0f}s]", flush=True)


def edges_of(pts, zs):
    cell = {}
    for i, z in enumerate(zs):
        cell.setdefault((int(z.real // 1), int(z.imag // 1)), []).append(i)
    out = []
    for i, z in enumerate(zs):
        cx, cy = int(z.real // 1), int(z.imag // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i or abs(abs(z - zs[j]) - 1) > 1e-6:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        out.append((i, j))
    return out


built = []
for ri, u in enumerate(ROT):
    assert knorm2(u) == K1
    img = [kmul(u, q) for q in Sa]
    imset = set(img)
    pts, seen2 = list(Sa), set(Sa)
    for q in img:
        if q not in seen2:
            seen2.add(q)
            pts.append(q)
    zs = [zof(q) for q in pts]
    E = edges_of(pts, zs)
    only_a, only_b = saset - imset, imset - saset
    cross = sum(1 for i, j in E
                if (pts[i] in only_a and pts[j] in only_b)
                or (pts[j] in only_a and pts[i] in only_b))
    built.append((cross, len(saset & imset), ri, pts, zs, E))
    if ri % 10 == 0:
        print(f"  ... built {ri}/{len(ROT)}  [{time.time()-t0:.0f}s]",
              flush=True)
built.sort(key=lambda t: (-t[0], t[1]))
print(f"\ncross-edge counts: {[b[0] for b in built[:20]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

for cross, shared, ri, pts, zs, E in built:
    if cross == 0:
        print(f"  the rest have no cross edge; stopping  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        break
    n = len(pts)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** rotation {ri}: {n} pts {len(E)} edges NOT 4-COLOURABLE "
              f"*** [{time.time()-t0:.0f}s]", flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/kchi5.pkl", "wb") as fh:
            pickle.dump(([flat(q) for q in pts], E), fh)
        break
    hits, hard, seenc = [], [], 0
    for i in range(n):
        for j in range(i + 1, n):
            vv = abs(zs[i] - zs[j]) ** 2
            if vv > 40:
                continue
            D = Fr(round(vv * 1584), 1584)
            if abs(float(D) - vv) > 1e-6 or D == 1 \
                    or not closable_over(D, KCL):
                continue
            if knorm2(ksub(pts[j], pts[i])) != krat(D):
                continue
            seenc += 1
            sv.conf_budget(20000)
            r = sv.solve_limited(assumptions=[1 + i * 4, -(1 + j * 4)])
            if r is False:
                hits.append((i, j, D))
                print(f"  *** rotation {ri}: FORCED PAIR {i},{j} at D = {D} "
                      f"*** [{time.time()-t0:.0f}s]", flush=True)
            elif r is None:
                hard.append((i, j, D))
    for i, j, D in hard:
        if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
            hits.append((i, j, D))
            print(f"  *** rotation {ri}: FORCED (full) {i},{j} at D = {D} ***",
                  flush=True)
    sv.delete()
    print(f"  rotation {ri}: {n} pts, {len(E)} edges, {cross} cross, "
          f"{shared} shared, {seenc} pairs, {len(hits)} forced, {len(hard)} "
          f"hard  [{time.time()-t0:.0f}s]", flush=True)
    if hits:
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/kforced.pkl", "wb") as fh:
            pickle.dump((ri, [flat(q) for q in pts], E, hits), fh)
        break
