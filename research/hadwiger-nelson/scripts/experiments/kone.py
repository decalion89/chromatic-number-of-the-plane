"""One biting rotation, verified exactly, then scanned.

The sweep counts cross edges by float, and the elements of K it works on carry
long coordinates, so the count has to be confirmed in exact arithmetic before
anything is read off it.  Take a handful of the rotations the sweep flagged,
build the union with exact unit-distance tests, and report what is really
there: how many points the copies share, how many edges cross, and whether
any pair at a K-closable distance is forced.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
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
saset = set(Sa)
print(f"Sa over K: {len(Sa)} points  [{time.time()-t0:.0f}s]", flush=True)


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


for ri in (372, 542, 721, 887, 1014, 1278):
    u = steps[ri]
    img = [kmul(u, q) for q in Sa]
    shared = len(set(img) & saset)
    pts, seen2 = list(Sa), set(Sa)
    for q in img:
        if q not in seen2:
            seen2.add(q)
            pts.append(q)
    zs = [zof(q) for q in pts]
    E = edges_of(pts, zs)
    n = len(pts)
    only_a = saset - set(img)
    only_b = set(img) - saset
    cross = sum(1 for i, j in E
                if (pts[i] in only_a and pts[j] in only_b)
                or (pts[j] in only_a and pts[i] in only_b))
    print(f"\nrotation {ri}: {n} points, {shared} shared, {len(E)} edges, "
          f"{cross} EXACT cross edges  [{time.time()-t0:.0f}s]", flush=True)
    if cross == 0:
        print("  the float count was spurious -- nothing to scan", flush=True)
        continue
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** NOT 4-COLOURABLE ***", flush=True)
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
                print(f"  *** FORCED PAIR {i},{j} at D = {D} ***", flush=True)
            elif r is None:
                hard.append((i, j, D))
    for i, j, D in hard:
        if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
            hits.append((i, j, D))
            print(f"  *** FORCED (full) {i},{j} at D = {D} ***", flush=True)
    sv.delete()
    print(f"  {seenc} pairs scanned, {len(hits)} forced, {len(hard)} hard  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if hits:
        break
