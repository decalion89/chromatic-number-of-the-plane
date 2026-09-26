"""Every ring that bites, and every pair it pins -- not just the antipodal one.

Two things loosen the search at once.

The antipodal pair is what de Grey's construction aims at, but nothing says a
union has to force THAT pair.  Any forced pair at a closable distance can be
spindled, so the second radical, sqrt(16D-1), is a condition on whatever pair
turns up, not on the ring chosen beforehand.  Dropping it from the ring test
takes G*'s candidates from three to twelve.

And a rotation has to bite.  rho_16 produces zero cross edges on G* -- it maps
Ya onto Yb, being the very angle de Grey's two turns differ by, so it extends
his fan rather than opening anything.  A union with no cross edges cannot
force what its halves do not already force, and G* forces nothing.

So: every rational closable ring, ranked by the cross edges its rotation
actually makes, and on each union the full agreeing-pair filter -- which
bucketing makes free, so all 27000-odd points are interrogated rather than a
handful of ring pairs.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import Counter
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point, _rot60, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance, agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/hn/")
ONE = K.rational(1)
t0 = time.time()
PIV = Point(K.rational(-2), K.zero())
rot60 = _rot60(K).about(PIV)

G = build_G(K, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot60(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
GE = sorted((min(a, b), max(a, b)) for a, b in build_graph(Gs).edges())
print(f"G*: {len(Gs)} points, {len(GE)} edges  [{time.time()-t0:.0f}s]",
      flush=True)

rings = Counter()
for p in Gs:
    dx, dy = p.x - PIV.x, p.y - PIV.y
    rings[float(dx * dx + dy * dy)] += 1
cand = []
for r2f, cnt in rings.items():
    D = Fr(round(r2f * 55440), 55440)
    if abs(float(D) - r2f) > 1e-9 or D < Fr(1, 4) or D == 1:
        continue
    if closable_distance(D):
        cand.append((cnt, D))
cand.sort(reverse=True)
print(f"{len(cand)} rational closable rings: "
      f"{[(c, str(d)) for c, d in cand]}  [{time.time()-t0:.0f}s]",
      flush=True)


def turn(pts, E, rot):
    out, idx = list(pts), {p: i for i, p in enumerate(pts)}
    img = []
    for p in pts:
        q = rot(p)
        j = idx.get(q)
        if j is None:
            j = len(out)
            idx[q] = j
            out.append(q)
        img.append(j)
    ed = set(E)
    for a, b in E:
        x, y = img[a], img[b]
        ed.add((min(x, y), max(x, y)))
    zf = [(float(p.x), float(p.y)) for p in out]
    cell = {}
    for i, (a, b) in enumerate(zf):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    cross = 0
    for j in range(len(pts), len(out)):
        a, b = zf[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if i >= j or abs((a - zf[i][0]) ** 2
                                     + (b - zf[i][1]) ** 2 - 1) > 1e-7:
                        continue
                    dx, dy = out[i].x - out[j].x, out[i].y - out[j].y
                    if dx * dx + dy * dy == ONE and (i, j) not in ed:
                        ed.add((i, j))
                        cross += 1
    return out, sorted(ed), cross, zf


results = []
for cnt, D in cand:
    U, UE, cross, zf = turn(Gs, GE, rotation_joining(D, K).about(PIV))
    n = len(U)
    print(f"  D={D} ({cnt} on the ring): {n} points, {len(UE)} edges, "
          f"{cross} CROSS  [{time.time()-t0:.0f}s]", flush=True)
    results.append((cross, D, n))
    if cross == 0:
        print("      no cross edges -- the rotation is internal, skipping",
              flush=True)
        continue
    k = 5
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in UE:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** D={D}: NOT 5-COLOURABLE ***", flush=True)
        with open(SC + "ringscan_six.pkl", "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in U], UE), fh)
        sv.delete()
        break
    rng, cols = random.Random(16180), []
    for s in range(12):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + v)
                       for v in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    pairs = agreeing_pairs(cols, colours=k)
    good = []
    for i, j in pairs:
        v = (zf[i][0] - zf[j][0]) ** 2 + (zf[i][1] - zf[j][1]) ** 2
        Dp = Fr(round(v * 55440), 55440)
        if abs(float(Dp) - v) > 1e-8 or Dp == 1 or Dp < Fr(1, 4):
            continue
        if closable_distance(Dp):
            good.append((i, j, Dp))
    forced = []
    for i, j, Dp in good:
        if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)]):
            forced.append((i, j, str(Dp)))
    sv.delete()
    print(f"      {len(pairs)} pairs agree in all 12 samples, {len(good)} at "
          f"a closable distance, {len(forced)} FORCED  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if forced:
        with open(SC + f"ringscan_forced_{D.numerator}_{D.denominator}.pkl",
                  "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in U], UE,
                         forced), fh)
        print(f"  *** FORCED PAIR AT FIVE COLOURS: {forced[:3]} ***",
              flush=True)
        break
print(f"cross by ring: {sorted(results, reverse=True)}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
