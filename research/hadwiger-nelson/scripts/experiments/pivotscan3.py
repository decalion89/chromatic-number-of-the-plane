"""Diversify by ring, not by pivot -- the ranking was returning one rotation.

Ranking 12199 (pivot, ring) unions by cross edges puts the same ring on top
twenty times over: D = 1/3, whose rotation has cos = 1 - 1/(2/3) = -1/2 and is
the 120-degree turn.  Different pivots, one rotation.  Testing the top sixty
would test that turn sixty times and every other closable rotation zero times,
which is not a scan of the space but a scan of one point in it.

So take the best union for each DISTINCT ring, and then the next best, and so
on -- breadth across rotations first, depth within a rotation only after.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import defaultdict
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance, agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/hn/")
ONE = K.rational(1)
PER_RING = 3
t0 = time.time()
g = build_G(K)
P = list(g.vertices)
GE = sorted((min(a, b), max(a, b)) for a, b in g.edges())
zf = [(float(p.x), float(p.y)) for p in P]
with open(SC + "pivotscan2.pkl", "rb") as fh:
    scored = pickle.load(fh)
by_ring = defaultdict(list)
for cross, v, Ds, _ in scored:
    if cross > 0:
        by_ring[Ds].append((cross, v))
order = []
for Ds, lst in by_ring.items():
    lst.sort(reverse=True)
    for cross, v in lst[:PER_RING]:
        order.append((cross, v, Fr(Ds)))
order.sort(reverse=True)
print(f"{len(by_ring)} distinct rings bite; {len(order)} unions to filter "
      f"(up to {PER_RING} per ring).  Best per ring: "
      f"{sorted(((max(l)[0], d) for d, l in by_ring.items()), reverse=True)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def turn(rot):
    out, idx = list(P), {p: i for i, p in enumerate(P)}
    zz, img = list(zf), []
    for p in P:
        q = rot(p)
        j = idx.get(q)
        if j is None:
            j = len(out)
            idx[q] = j
            out.append(q)
            zz.append((float(q.x), float(q.y)))
        img.append(j)
    ed = set(GE)
    for a, b in GE:
        x, y = img[a], img[b]
        ed.add((min(x, y), max(x, y)))
    cell = {}
    for i, (a, b) in enumerate(zz):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    cross = 0
    for j in range(len(P), len(out)):
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


k = 5
best = (-1, None)
for nth, (approx, v, D) in enumerate(order):
    U, UE, cross, zz = turn(rotation_joining(D, K).about(P[v]))
    n = len(U)
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in UE:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** pivot {v} D={D}: {n} points, {cross} cross -- NOT "
              f"5-COLOURABLE ***", flush=True)
        with open(SC + "pivotscan3_six.pkl", "wb") as fh:
            pickle.dump((v, str(D), [(str(p.x), str(p.y)) for p in U], UE),
                        fh)
        sv.delete()
        break
    rng, cols = random.Random(99991 + nth), []
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
    if len(pairs) > best[0]:
        best = (len(pairs), f"pivot {v} D={D}")
    print(f"  D={D} pivot {v}: {n} points, {cross} cross, {len(pairs)} agree,"
          f" {len(good)} closable, {len(forced)} FORCED   (best agreement so "
          f"far: {best[0]} at {best[1]})  [{time.time()-t0:.0f}s]",
          flush=True)
    if forced:
        with open(SC + "pivotscan3_forced.pkl", "wb") as fh:
            pickle.dump((v, str(D), [(str(p.x), str(p.y)) for p in U], UE,
                         forced), fh)
        print(f"  *** FORCED PAIR AT FIVE COLOURS {forced[:3]} ***",
              flush=True)
        break
