"""Drive down the number of critical vertices, which lower-bounds rho.

rho, the least set using all k colours in every colouring, is at least the
number of CRITICAL vertices: for each u with W - u still (k-1)-colourable,
that colouring leaves u alone in the kth colour, so every forcing set contains
u. In de Grey's G all 1581 qualify, which is exactly why rho = n there.

Adding points can only help. A vertex stops being critical as soon as W - u
needs k colours anyway, and extra points only make colouring harder. So the
count is a monotone gradient -- the first one in this package pointing
directly at rho -- and the natural points to add are the graph's own holes,
the exact points one away from two or more vertices at once.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.graph import build_graph
from hn.mixed import circle_intersections

g = build_G()
one = g.vertices[0].x.field.rational(1)
xs = [float(v.x) for v in g.vertices]
ys = [float(v.y) for v in g.vertices]
cells = collections.defaultdict(list)
for i, (a, b) in enumerate(zip(xs, ys)):
    cells[(int(a // 2), int(b // 2))].append(i)

t0 = time.time()
have = set(g.vertices)
holes = collections.Counter()
for i in range(g.n):
    cx, cy = int(xs[i] // 2), int(ys[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i or (xs[i]-xs[j])**2 + (ys[i]-ys[j])**2 > 3.999:
                    continue
                for q in circle_intersections(g.vertices[i], one,
                                              g.vertices[j], one):
                    if q not in have:
                        holes[q] += 1
    if i % 500 == 499:
        print(f"    ... {i+1}/{g.n}, {len(holes)} holes  "
              f"[{time.time()-t0:.0f}s]", flush=True)
ranked = sorted(holes, key=lambda q: -holes[q])
print(f"{len(holes)} exact holes; the best is one away from "
      f"{holes[ranked[0]]} vertices  [{time.time()-t0:.0f}s]", flush=True)

random.seed(43)
for take in (0, 500, 2000, 6000, len(ranked)):
    pts = list(g.vertices) + ranked[:take]
    w = build_graph(pts)
    if not is_k_colorable(w, 5)[0]:
        print(f"  *** +{take} holes: {w} is NOT 5-COLOURABLE -- "
              f"chi(R^2) >= 6 ***", flush=True)
        break
    sample = random.sample(range(w.n), 12)
    crit = 0
    for u in sample:
        sub = build_graph([q for m, q in enumerate(w.vertices) if m != u])
        if is_k_colorable(sub, 4)[0]:
            crit += 1
    print(f"  +{take} holes: {w}  critical on {crit}/{len(sample)} sampled "
          f"vertices  [{time.time()-t0:.0f}s]", flush=True)
