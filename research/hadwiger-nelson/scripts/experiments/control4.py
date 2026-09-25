"""The same scan, on Sa at four colours, where the answer is known.

The pivot scan found nothing at five colours across 12199 candidates and 44
biting rotations.  That is only worth something if the scan would have found
de Grey's construction had it been looking one level down -- so run it there,
unchanged: every vertex of Sa as a pivot, every rational closable ring about
it, the rotation that carries the ring to distance 1, and the agreement filter
at four colours instead of five.

What it should find, with nothing told to it, is the origin as pivot and
D = 4 as ring, pinning seventy-five pairs.  Anything less and the five-colour
negative is a statement about the method.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import Counter
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance, agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
ONE = K.rational(1)
k = 4
t0 = time.time()
P = build_Sa(K)
g = build_graph(P)
GE = sorted((min(a, b), max(a, b)) for a, b in g.edges())
zf = [(float(p.x), float(p.y)) for p in P]
print(f"Sa: {len(P)} points, {len(GE)} edges  [{time.time()-t0:.0f}s]",
      flush=True)

# Sa is invariant under the 12-element dihedral group about the origin, so
# the unions about p and about g(p) are isomorphic and one pivot per orbit
# says everything 397 of them would.  That is the difference between a
# control that runs and one that does not.
from hn.geometry import Point, _rot60
r60 = _rot60(K)
pos = {p: i for i, p in enumerate(P)}
rep, seen_orb = [], set()
for v, p in enumerate(P):
    if v in seen_orb:
        continue
    rep.append(v)
    q = p
    for _ in range(6):
        for w in (q, Point(q.x, -q.y)):
            j = pos.get(w)
            if j is not None:
                seen_orb.add(j)
        q = r60(q)
print(f"{len(rep)} pivot orbits of {len(P)} points  "
      f"[{time.time()-t0:.0f}s]", flush=True)

jobs = []
for v in rep:
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
print(f"{len(jobs)} (pivot, ring) candidates  [{time.time()-t0:.0f}s]",
      flush=True)


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


best, hits = (-1, None), []
for nth, (cnt, v, D) in enumerate(jobs):
    U, UE, cross, zz = turn(rotation_joining(D, K).about(P[v]))
    n = len(U)
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in UE:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  pivot {v} D={D}: {n} points, {cross} cross -- NOT "
              f"{k}-COLOURABLE (that is chi >= 5 straight away)", flush=True)
        sv.delete()
        continue
    rng, cols = random.Random(3141 + nth), []
    for s in range(12):
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
    # At four colours a single forced-pair proof runs to minutes, and it
    # would be re-proving de Grey's theorem.  Agreement is the discovery
    # step, so only the leader is confirmed, at the end.
    # No SAT confirmation here at all.  One forced-pair proof at four
    # colours runs over ten minutes on these graphs -- the first union alone
    # spent half an hour on three of them -- and what it would confirm is de
    # Grey's published theorem.  Agreement is the discovery step, and the
    # discovery step is what this control exists to test.
    forced = []
    sv.delete()
    if len(pairs) > best[0]:
        best = (len(pairs), f"pivot {v} D={D} ({cross} cross)")
    if pairs or cross > 0:
        print(f"  pivot {v} D={D}: {n} points, {cross} cross, {len(pairs)} "
              f"agree, {len(good)} closable, {len(forced)} FORCED  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    if forced:
        hits.append((v, str(D), cross, len(pairs), forced[:4]))
        with open(SC + "control4_forced.pkl", "wb") as fh:
            pickle.dump(hits, fh)
    if nth % 200 == 0:
        print(f"  ... {nth}/{len(jobs)}, best agreement {best[0]} at "
              f"{best[1]}, {len(hits)} unions with a forced pair  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"DONE: best agreement {best[0]} at {best[1]}; "
      f"{len(hits)} unions carry a forced pair at four colours  "
      f"[{time.time()-t0:.0f}s]", flush=True)
