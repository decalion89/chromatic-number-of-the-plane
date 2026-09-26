"""Where must the 60-degree rotation be centred for G's orbit to overlap?

de Grey's Sa is the dihedral group D6 applied to S about the origin -- the
60-degree rotation and the reflection in y, twelve isometries -- and it is
397 points rather than 39*12 = 468, so the orbit folds onto itself 71 times.
That folding is where its rigidity comes from, and it is why Sa reaches
pressure 3 at four colours while a single translated copy of G reaches only
2 at five.

Applying the same group to G about the origin gives 18966 points, exactly
1581*12: no coincidences at all, so the twelve copies are disjoint and
constrain each other only through whatever edges happen to fall between them.
The centre is wrong, not the construction.

The right centres are computable. A rotation by theta about c takes u to v
exactly when c = (v - R u)/(I - R), and for theta = 60 degrees the matrix
I - R has determinant 1 with inverse R, so c = R v - R^2 u -- one centre per
ordered pair, and the centres that recur most are the ones whose orbit folds.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from math import cos, sin, pi
from hn.degrey import build_G

g = build_G()
xs = [float(v.x) for v in g.vertices]
ys = [float(v.y) for v in g.vertices]
c60, s60 = cos(pi / 3), sin(pi / 3)


def R(x, y):
    return c60 * x - s60 * y, s60 * x + c60 * y


t0 = time.time()
cnt, rep = collections.Counter(), {}
SC = 10 ** 6
for i in range(g.n):
    ux, uy = R(*R(xs[i], ys[i]))          # R^2 u
    for j in range(g.n):
        vx, vy = R(xs[j], ys[j])          # R v
        cx, cy = vx - ux, vy - uy
        if cx * cx + cy * cy > 400.0:
            continue
        key = (round(cx * SC), round(cy * SC))
        cnt[key] += 1
        if key not in rep:
            rep[key] = (i, j)
    if i % 300 == 299:
        print(f"    ... {i+1}/{g.n}, {len(cnt)} centres  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"{len(cnt)} candidate centres  [{time.time()-t0:.0f}s]", flush=True)
for key, k in cnt.most_common(12):
    print(f"   centre ({key[0]/SC:+.5f}, {key[1]/SC:+.5f})  folds {k} points",
          flush=True)
pickle.dump([(k, key, rep[key]) for key, k in cnt.most_common(30)],
            open(f"{SC if False else '/tmp/hn'}/rotcentres.pkl", "wb"))
print("saved top 30", flush=True)
