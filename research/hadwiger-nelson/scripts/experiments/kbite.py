"""The pair scan, on the rotations that actually bite.

de Grey's Y is Sa and its image sharing one point and joined by six edges, and
that is enough to turn a graph with no forced pair into one that has one.  The
sweep finds rotations of K doing the same thing harder -- one shared point and
thirty cross edges -- so those are where the scan belongs.

A forced pair at a K-closable distance here is a 5-chromatic unit-distance
graph over the field that blocks at the gate, which no multiquadratic field
can carry.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/hn/"
           "kneck2.py").read()
exec(src[:src.index("pairsum = {}")])
from hn.degrey import S_POINTS
from hn.homcol import closable_over
from pysat.solvers import Solver

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
        q = base
        for _ in range(6):
            if q not in seenp:
                seenp.add(q)
                Sa.append(q)
            q = kmul(Z6, q)
with open("/tmp/hn/kcross.pkl", "rb") as fh:
    raw = pickle.load(fh)


def unflat(t):
    return tuple(tuple(t[3 * i:3 * i + 3]) for i in range(4))


rot = [(c, s, unflat(v)) for c, s, v in raw]
rot.sort(key=lambda r: (-r[0], r[1]))
print(f"{len(rot)} biting rotations, best {rot[0][0]} cross edges with "
      f"{rot[0][1]} shared  [{time.time()-t0:.0f}s]", flush=True)


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
                    if j <= i or abs(abs(z - zs[j]) - 1) > 1e-7:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        out.append((i, j))
    return out


seen_rot = set()
for ri, (cross, shared, u) in enumerate(rot):
    key = min(tuple(flat(kmul(u, kmul(Z6, K1) if False else K1))) for _ in [0])
    orb = []
    v = u
    for _ in range(6):
        orb.append(flat(v))
        v = kmul(v, Z6)
    key = min(orb)
    if key in seen_rot:
        continue
    seen_rot.add(key)
    pts, seen2 = list(Sa), set(Sa)
    for q in Sa:
        z = kmul(u, q)
        if z not in seen2:
            seen2.add(z)
            pts.append(z)
    zs = [zof(q) for q in pts]
    E = edges_of(pts, zs)
    n = len(pts)
    cls = [[1 + v2 * 4 + c for c in range(4)] for v2 in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** rotation {ri}: {n} pts {len(E)} edges NOT 4-COLOURABLE "
              f"*** [{time.time()-t0:.0f}s]", flush=True)
        with open("/tmp/hn/kchi5.pkl", "wb") as fh:
            pickle.dump(([flat(q) for q in pts], E), fh)
        break
    hits, hard, seenc = [], [], 0
    for i in range(n):
        for j in range(i + 1, n):
            vv = abs(zs[i] - zs[j]) ** 2
            if vv > 40:
                continue
            D = Fr(round(vv * 1584), 1584)
            if abs(float(D) - vv) > 1e-7 or D == 1 \
                    or not closable_over(D, KCL):
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
    print(f"  rotation {ri} ({cross} cross, {shared} shared): {n} pts, "
          f"{len(E)} edges, {seenc} pairs, {len(hits)} forced, {len(hard)} "
          f"hard  [{time.time()-t0:.0f}s]", flush=True)
    if hits:
        with open("/tmp/hn/kforced.pkl", "wb") as fh:
            pickle.dump((ri, [flat(q) for q in pts], E, hits), fh)
        break
