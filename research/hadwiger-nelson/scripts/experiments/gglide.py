"""Glide reflections: the last family of plane isometries left untried.

Rotations were solved for, translates were histogrammed, and both move the
graph rigidly forward.  A reflection does not: sigma(z) = u.zbar reverses
orientation, so sigma(G) is a genuinely different point set from every
rotation and every translate of G -- and G is not symmetric under it, Y not
being.

The equation is unchanged.  A cross edge for sigma followed by a translate is
|p - (sigma(q) + t)| = 1, which is |p - qbar' - t| = 1 with qbar' = sigma(q):
the same histogram of t = p - sigma(q) - v, over the conjugated point set.
So the machinery ports with one conjugation.


Every union here has been G with a ROTATED copy.  But G u (G + t) is just as
good a unit-distance graph, just as 5-chromatic, and translates cost nothing
arithmetically: a cross edge is |p - (q + t)| = 1, i.e.

    t = p - q - v      for p, q in G and v a unit direction of the plane,

so every difference of G, shifted by any unit vector, is a candidate.  No
square has to lie in the field, no rotation has to exist.  The family is far
larger than the rotations and it has gone completely unexamined.

What makes a translate good is multiplicity: the number of (p, q, v) with
p - q - v = t is exactly the number of cross edges it buys.  Differences of G
concentrate heavily -- it is built from rotated copies of a small seed -- so
the histogram is worth computing before anything else.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from collections import Counter
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
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
P = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c)) for p in pts]
zf = [(float(p.x), float(p.y)) for p in pts]
GE = list(build_graph(pts).edges())
ONE = (Fr(1),) + (Fr(0),) * (DIM - 1)


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


DIRS = set()
for a, b in GE:
    d = (msub(P[b][0], P[a][0]), msub(P[b][1], P[a][1]))
    DIRS.add(d)
    DIRS.add((tuple(-c for c in d[0]), tuple(-c for c in d[1])))
DIRS = sorted(DIRS)
RAD = []
for m in range(DIM):
    pr = 1
    for bb in range(4):
        if m >> bb & 1:
            pr *= GENS[bb]
    RAD.append(pr ** 0.5)
print(f"G: {len(P)} points, {len(GE)} edges, {len(DIRS)} unit directions  "
      [0] if False else f"G: {len(P)} points, {len(GE)} edges, "
      f"{len(DIRS)} unit directions  [{time.time()-t0:.0f}s]", flush=True)

# The exact histogram costs two gigabytes in two minutes -- a translate is 32
# rationals and there are hundreds of thousands of them per point.  Key on the
# float pair instead, rounded well inside the separation of distinct elements
# of this field, and confirm the survivors exactly.  The coordinates here are
# short (small rationals over sqrt3, sqrt5, sqrt7, sqrt11), so the float is a
# safe filter -- unlike over K, where eighty-digit coordinates made a 1e-9
# tolerance report cross edges that did not exist.
DIRZ = []
for v in DIRS:
    DIRZ.append((sum(float(x) * r for x, r in zip(v[0], RAD)),
                 sum(float(x) * r for x, r in zip(v[1], RAD))))
step = max(1, len(P) // 150)
hist = Counter()
# sigma is reflection in the x-axis: (x, y) -> (x, -y).  Composed with every
# translate this sweeps the whole orientation-reversing half of the isometry
# group, which no search here has touched.
zs = [(x, -y) for x, y in zf]
for pi in range(0, len(P), step):
    px, py = zf[pi]
    for vx, vy in DIRZ:
        ax, ay = px - vx, py - vy
        for qx, qy in zs:
            hist[(round(ax - qx, 7), round(ay - qy, 7))] += 1
    if (pi // step) % 30 == 0:
        print(f"  ... p {pi}/{len(P)}, {len(hist)} translates, best "
              f"{hist.most_common(1)[0][1] if hist else 0}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
cands = [t for t, c in hist.most_common(400) if abs(t[0]) + abs(t[1]) > 1e-9]
print(f"\n{len(hist)} distinct translates; top counts "
      f"{[c for _, c in hist.most_common(12)]}  [{time.time()-t0:.0f}s]",
      flush=True)

# now the exact translate for each surviving float key: t = p - v - q for the
# (p, v, q) that produced it, recovered by searching the sample again
want = set(cands)
exact = {}
for pi in range(0, len(P), step):
    px, py = zf[pi]
    p = P[pi]
    for vi, (vx, vy) in enumerate(DIRZ):
        ax, ay = px - vx, py - vy
        v = DIRS[vi]
        for qi, (qx, qy) in enumerate(zs):
            k = (round(ax - qx, 7), round(ay - qy, 7))
            if k in want and k not in exact:
                q = (P[qi][0], tuple(-c for c in P[qi][1]))
                exact[k] = (msub(msub(p[0], v[0]), q[0]),
                            msub(msub(p[1], v[1]), q[1]))
    if len(exact) == len(want):
        break
print(f"{len(exact)} of {len(want)} recovered exactly  "
      f"[{time.time()-t0:.0f}s]", flush=True)
out = [(hist[k], exact[k]) for k in cands if k in exact]
out.sort(key=lambda t: -t[0])
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/gglide.pkl", "wb") as fh:
    pickle.dump(out, fh)
print(f"saved {len(out)} translates, best sampled count {out[0][0]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
