"""Back down to G's own size, and scan exhaustively instead of deeply.

Y has 793 points.  Every object built in this session to beat it has been
larger -- 3000, 5900, 13873, 27673, 55345 -- and none forces anything, while
the sampling that decides the question costs over eleven minutes a colouring
at 41509 points.  Size has been bought at the price of throughput, and size
was never the variable.

So come back to G's scale, where a colouring costs milliseconds, and spend the
budget on breadth.  The ring criterion needs no global symmetry: about ANY
pivot p, a ring of squared radius D carries its points to distance 1 from
their images under rotation_joining(D), and D only has to be rational with
sqrt(4D-1) in the field.  Every vertex of G is a candidate pivot and each has
its own rings, which is thousands of unions of about three thousand points
each rather than one of thirty thousand.

Cross edges are counted first, for everything, because counting is cheap and
solving is not; the filter then runs on the unions that actually bite.
"""
import sys, time, pickle, random
from collections import Counter
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance, agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
ONE = K.rational(1)
t0 = time.time()
g = build_G(K)
P = list(g.vertices)
GE = sorted((min(a, b), max(a, b)) for a, b in g.edges())
deg = [len(a) for a in g.adj]
zf = [(float(p.x), float(p.y)) for p in P]
print(f"G: {len(P)} points, {len(GE)} edges, max degree {max(deg)}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

order = sorted(range(len(P)), key=lambda v: -deg[v])   # every pivot


def turn(pts, E, zfl, rot):
    out, idx = list(pts), {p: i for i, p in enumerate(pts)}
    zz = list(zfl)
    img = []
    for p in pts:
        q = rot(p)
        j = idx.get(q)
        if j is None:
            j = len(out)
            idx[q] = j
            out.append(q)
            zz.append((float(q.x), float(q.y)))
        img.append(j)
    ed = set(E)
    for a, b in E:
        x, y = img[a], img[b]
        ed.add((min(x, y), max(x, y)))
    cell = {}
    for i, (a, b) in enumerate(zz):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    cross = 0
    for j in range(len(pts), len(out)):
        a, b = zz[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if i >= j or abs((a - zz[i][0]) ** 2
                                     + (b - zz[i][1]) ** 2 - 1) > 1e-7:
                        continue
                    dx, dy = out[i].x - out[j].x, out[i].y - out[j].y
                    if dx * dx + dy * dy == ONE and (i, j) not in ed:
                        ed.add((i, j))
                        cross += 1
    return out, sorted(ed), cross, zz


# Phase one: count, do not solve.
jobs = []
for nth, v in enumerate(order):
    px, py = zf[v]
    rings = Counter()
    for i, (a, b) in enumerate(zf):
        if i != v:
            rings[round((a - px) ** 2 + (b - py) ** 2, 9)] += 1
    for r2f, cnt in rings.items():
        D = Fr(round(r2f * 55440), 55440)
        if abs(float(D) - r2f) > 1e-8 or D < Fr(1, 4) or D == 1:
            continue
        if closable_distance(D):
            jobs.append((cnt, v, D))
    if nth % 60 == 0:
        print(f"  ... pivots {nth}/{len(order)}, {len(jobs)} (pivot, ring) "
              f"pairs  [{time.time()-t0:.0f}s]", flush=True)
jobs.sort(reverse=True)
print(f"{len(jobs)} (pivot, ring) candidates; top populations "
      f"{[c for c, _, _ in jobs[:12]]}  [{time.time()-t0:.0f}s]", flush=True)

# Scoring in floating point.  Counting is a ranking, not a verdict: the
# coordinates are O(10) with small denominators, so double precision carries
# about 1e-14 against a 1e-7 window -- seven orders of margin, the same
# argument build_graph itself rests on -- and every union that reaches phase
# two is rebuilt exactly anyway.  Rotating 1581 points through a
# 16-dimensional field costs 1.6 million rational operations and a second and
# a half; in floats it is microseconds, which is the difference between 240
# pivots and all of them.
# Scoring, vectorised.  The inner loop over 1581 points was the cost, not
# the grid: Python pays per point, numpy pays per array.  The whole distance
# matrix between a rotated copy and the original is 1581 x 1581 -- 2.5 million
# entries, twenty megabytes, twenty milliseconds -- so the grid is not worth
# keeping at all at this size.  Twelve thousand unions in five minutes instead
# of eighty-five.
import numpy as np
XY = np.array(zf)
RJ = {}

scored = []
for nth, (cnt, v, D) in enumerate(jobs):
    if D not in RJ:
        r = rotation_joining(D, K)
        RJ[D] = (float(r.cos), float(r.sin))
    c, sn = RJ[D]
    px, py = zf[v]
    dx, dy = XY[:, 0] - px, XY[:, 1] - py
    qx, qy = px + c * dx - sn * dy, py + sn * dx + c * dy
    d2 = ((qx[:, None] - XY[None, :, 0]) ** 2
          + (qy[:, None] - XY[None, :, 1]) ** 2)
    # An image that lands on an existing point contributes G's own edges, not
    # cross edges.  Dropping those rows is what keeps a near-symmetry -- a
    # rotation like rho_16, which maps half the graph onto the other half --
    # from scoring highest while adding nothing.
    fresh = d2.min(axis=1) > 1e-16
    cross = int(((np.abs(d2 - 1.0) < 1e-7) & fresh[:, None]).sum())
    scored.append((cross, v, D, 0))
    if nth % 1000 == 0:
        print(f"  ... scoring {nth}/{len(jobs)}, best {max(scored)[0]} cross "
              f"[{time.time()-t0:.0f}s]", flush=True)

scored.sort(reverse=True)
with open(SC + "pivotscan2.pkl", "wb") as fh:
    pickle.dump([(c, v, str(D), n) for c, v, D, n in scored], fh)
print(f"cross edges, top 20: {[(c, str(D)) for c, _, D, _ in scored[:20]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

# Phase two: solve the ones that bite hardest.
k = 5
for nth, (cross, v, D, _) in enumerate(scored[:60]):
    if cross == 0:
        break
    U, UE, _, zz = turn(P, GE, zf, rotation_joining(D, K).about(P[v]))
    n = len(U)
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in UE:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** pivot {v} D={D}: {n} points, {cross} cross -- NOT "
              f"5-COLOURABLE ***", flush=True)
        with open(SC + "pivotscan2_six.pkl", "wb") as fh:
            pickle.dump((v, str(D), [(str(p.x), str(p.y)) for p in U], UE),
                        fh)
        sv.delete()
        break
    rng, cols = random.Random(57721 + nth), []
    for s in range(14):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    pairs = agreeing_pairs(cols, colours=k)
    good = []
    for i, j in pairs:
        val = (zz[i][0] - zz[j][0]) ** 2 + (zz[i][1] - zz[j][1]) ** 2
        Dp = Fr(round(val * 55440), 55440)
        if abs(float(Dp) - val) > 1e-8 or Dp == 1 or Dp < Fr(1, 4):
            continue
        if closable_distance(Dp):
            good.append((i, j, Dp))
    forced = [(i, j, str(Dp)) for i, j, Dp in good
              if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    sv.delete()
    print(f"  pivot {v} D={D}: {n} points, {len(UE)} edges, {cross} cross, "
          f"{len(pairs)} agree, {len(good)} closable, {len(forced)} FORCED  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if forced:
        with open(SC + "pivotscan2_forced.pkl", "wb") as fh:
            pickle.dump((v, str(D), [(str(p.x), str(p.y)) for p in U], UE,
                         forced), fh)
        print(f"  *** FORCED PAIR AT FIVE COLOURS {forced[:3]} ***",
              flush=True)
        break
