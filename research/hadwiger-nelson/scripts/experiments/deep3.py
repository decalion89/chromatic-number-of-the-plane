"""Go straight to the depth the law points at, and sample rather than sweep.

Fitting dearest = 28 . exp(cross/300) to the control (0 cross, 28), a glide
reflection (393, 151) and a translate (1442, 3689) puts the forced regime --
where Y's own forced pair sits, above two million conflicts -- at about 3350
cross edges.  Depth 2 of the translate stack carries 2864 and depth 3 should
carry near 4300.

So build depth 3 and 4 directly.  And sample the pairs instead of sweeping
them: if the law holds, forced pairs are not rare at that depth, and a few
hundred queries will meet one.  A sweep of eighty thousand at a hundred
thousand conflicts each is a day's work for the same answer.

Every query runs with no budget at all -- the budget lesson stands.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
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


with open("/tmp/hn/gtrans.pkl", "rb") as fh:
    TR = pickle.load(fh)
random.seed(20260920)

for depth in (3, 4):
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
    cross = len(E) - inside
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    five = sv.solve()
    print(f"depth {depth}: {n} points, {len(E)} edges, {cross} cross, "
          f"5-colourable {five}; the law predicts a dearest pair near "
          f"{int(28 * 2.718281828 ** (cross / 300)):.0f}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not five:
        print("  *** SIX COLOURS ***", flush=True)
        with open("/tmp/hn/deep_six.pkl", "wb") as fh:
            pickle.dump((depth, allp, E), fh)
        break
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
    random.shuffle(cand)
    sample = cand[:400]
    print(f"  {len(cand)} pairs at a closable distance; sampling "
          f"{len(sample)}  [{time.time()-t0:.0f}s]", flush=True)

    def conf():
        return sv.accum_stats().get("conflicts", 0)

    base, dear, hits = conf(), 0, []
    for k, (i, j) in enumerate(sample):
        r = sv.solve(assumptions=[1 + i * 5, -(1 + j * 5)])
        c2 = conf()
        dear = max(dear, c2 - base)
        base = c2
        if not r:
            hits.append((i, j))
            print(f"  *** FORCED PAIR AT FIVE COLOURS: {i},{j} ***  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            with open("/tmp/hn/deep_forced.pkl",
                      "wb") as fh:
                pickle.dump((depth, allp, E, (i, j)), fh)
            break
        if k and k % 50 == 0:
            print(f"    ... {k}/{len(sample)}, dearest so far {dear}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"  depth {depth}: {len(hits)} forced in {len(sample)} sampled, "
          f"dearest {dear}  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    if hits:
        break
