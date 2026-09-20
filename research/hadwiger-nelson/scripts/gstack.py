"""Stack copies about one pivot, which is what Sa is.

Sa is not one union: it is S under twelve rotations, and that is where its
tightness comes from.  A single pivot union of G buys 796 cross edges and does
not force, so the move is to stack -- G, u.G, u^2.G, u^3.G about the same
vertex, each power adding its own cross edges to every copy below it.

The powers of one rotation are free: u has infinite order unless it is a root
of unity, so u^k gives a genuinely new copy each time, and all of them share
the pivot.  Scan after each addition, because the point is to find where the
forcing starts, not to build the biggest graph.
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


paths = sorted(glob.glob("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-"
                         "5848-a506-39c59179b415/scratchpad/gpivot_*.pkl"))
with open(paths[0], "rb") as fh:
    pv, ROT = pickle.load(fh)
o = P0[pv]
P = [(msub(p[0], o[0]), msub(p[1], o[1])) for p in P0]
oz = (fl(o[0]), fl(o[1]))
zt = [(z[0] - oz[0], z[1] - oz[1]) for z in zf]
print(f"pivot {pv}, {len(ROT)} rotations; stacking powers  "
      f"[{time.time()-t0:.0f}s]", flush=True)

# ROT[0] turns out to be u = -1, the half turn about the pivot: u^2 = 1, so
# the second power adds nothing and the stack stalls at 2371 points.  Stack
# DISTINCT rotations instead, which is what Sa does -- twelve different
# elements of a group, not the powers of one.
for depth in (1, 2, 3, 4, 5, 6):
    allp, allz = list(P), list(zt)
    idx = {p: i for i, p in enumerate(P)}
    copies = [list(range(len(P)))]
    for k in range(depth):
        v = ROT[k % len(ROT)]
        vx, vy = fl(v[0]), fl(v[1])
        this = []
        for p, z in zip(P, zt):
            q = cmul(v, p)
            if q in idx:
                this.append(idx[q])
            else:
                idx[q] = len(allp)
                this.append(len(allp))
                allp.append(q)
                allz.append((vx * z[0] - vy * z[1], vy * z[0] + vx * z[1]))
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
    five = sv.solve()
    print(f"  depth {depth}: {n} points, {len(E)} edges "
          f"({len(E)-inside} cross), 5-colourable {five}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not five:
        print("  *** SIX COLOURS ***", flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/gsix.pkl", "wb") as fh:
            pickle.dump((pv, depth, allp, E), fh)
        break
    cand = []
    for i in range(n):
        ai, bi = allz[i]
        for j in range(i + 1, n):
            vv = (ai - allz[j][0]) ** 2 + (bi - allz[j][1]) ** 2
            if vv > 36.0:
                continue
            D = Fr(round(vv * 1584), 1584)
            if abs(float(D) - vv) > 1e-7 or D == 1 \
                    or not closable_distance(D):
                continue
            cand.append((i, j))
    # Conflicts, not just verdicts.  G alone costs 6410 over its whole
    # 21358-pair scan and its dearest pair 29; Y's forced pair costs over two
    # million.  If the rich unions cost the same as G, stacking is not
    # tightening anything, and that is worth knowing before stacking further.
    def conf():
        return sv.accum_stats().get("conflicts", 0)

    base = conf()
    dear = 0
    hits, hard = [], []
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
    sv.delete()
    print(f"    {len(cand)} pairs, {len(hits)} FORCED, {len(hard)} rerun, "
          f"{conf()} conflicts in total, dearest {dear} "
          f"(G alone: 6410 and 29)  [{time.time()-t0:.0f}s]", flush=True)
    if hits:
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/gforced.pkl", "wb") as fh:
            pickle.dump((pv, depth, allp, E, hits), fh)
        print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
              flush=True)
        break
