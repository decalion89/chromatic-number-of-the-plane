"""de Grey's move over K, about the point his own is about.

Sb is Sa turned about the ORIGIN, and Sa is dihedrally symmetric there, so
that union doubles a symmetry rather than gluing two lumps at a vertex.  The
scan so far has been over pivots ranked by ring size, which need not include
the origin at all -- and 35 unions at 9000 pairs each have turned up nothing.

So: every rotation K supplies, applied about the origin, and the pair scan at
every K-closable distance.  The rational-chord six are 1, 1/3, 3, 11/3, 11/9,
25/9, closing D = 1, 3, 1/3, 3/11, 9/11, 9/25 -- all K-closable, since
1 - 4D has squarefree part -3 or -11 in each case.  rho_7 = (13+3sqrt-3)/14
joins them, and a sample of the irrational-chord ones after that.

A forced pair here is a 5-chromatic unit-distance graph over the field that
blocks at the gate.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
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
zsA = [zof(p) for p in Sa]
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


ROTS = [(f"chord {c}", r) for c, r in rots]
ROTS.append(("rho_7", kmul(krat(Fr(1, 14)),
                           kadd(krat(13), kmul(krat(3), SQ3)))))
for i, u in enumerate(cands[:24]):
    ROTS.append((f"irrational {i}", u))
print(f"{len(ROTS)} rotations about the origin  [{time.time()-t0:.0f}s]",
      flush=True)

for name, rho in ROTS:
    if knorm2(rho) != K1:
        print(f"  {name}: not a unit, skipped", flush=True)
        continue
    pts, seen2 = list(Sa), set(Sa)
    for q in Sa:
        z = kmul(rho, q)
        if z not in seen2:
            seen2.add(z)
            pts.append(z)
    if len(pts) == len(Sa):
        print(f"  {name}: a symmetry of Sa, nothing new", flush=True)
        continue
    zs = [zof(q) for q in pts]
    E = edges_of(pts, zs)
    n = len(pts)
    cross = len(E) - 2 * 1974
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** {name}: {n} pts {len(E)} edges NOT 4-COLOURABLE ***",
              flush=True)
        sv.delete()
        break
    hits, hard, seenc = [], [], 0
    for i in range(n):
        for j in range(i + 1, n):
            v = abs(zs[i] - zs[j]) ** 2
            if v > 40:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_over(D, KCL):
                continue
            seenc += 1
            sv.conf_budget(20000)
            r = sv.solve_limited(assumptions=[1 + i * 4, -(1 + j * 4)])
            if r is False:
                hits.append((i, j, D))
                print(f"  *** {name}: FORCED PAIR {i},{j} at D = {D} ***",
                      flush=True)
            elif r is None:
                hard.append((i, j, D))
    for i, j, D in hard:
        if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
            hits.append((i, j, D))
            print(f"  *** {name}: FORCED (full) {i},{j} at D = {D} ***",
                  flush=True)
    sv.delete()
    print(f"  {name}: {n} pts, {len(E)} edges ({cross:+d} cross), {seenc} "
          f"pairs, {len(hits)} forced, {len(hard)} hard  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if hits:
        import pickle
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/kforced.pkl", "wb") as fh:
            pickle.dump((name, [flat(q) for q in pts], E, hits), fh)
        break
