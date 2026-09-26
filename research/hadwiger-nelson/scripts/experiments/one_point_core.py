"""G plus a single point: the sweet spot of the rigidity tension.

There is a real tension in what a core needs. Rigid graphs have few colourings
and so small cores -- but the most rigid are the vertex-critical ones, and
criticality kills cores outright, since a critical vertex is alone in the kth
colour in some colouring and so separable from everything at once. Loose
graphs escape criticality but their cores are huge: the unions let a
counterexample park the free colour anywhere, and 400 rounds of
counterexample-killing close nothing.

Adding ONE point to a 5-critical graph resolves it exactly. de Grey's G is
5-vertex-critical, so in W = G + v every OLD vertex u still has W - u
containing G - u, which is 4-colourable; only v itself has W - v = G, which is
5-chromatic. So v is the unique non-critical vertex of W and the unique pivot
that can carry a core -- while G keeps every bit of its rigidity.

The useful v are the exact points one away from several G-vertices: the more
neighbours, the more the circle constrains, and the pressure theorem then
bounds the core from below.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, itertools, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.forced import (ColourRelations, cegar_core, is_core, minimise_core,
                       pressure)
from hn.graph import build_graph
from hn.mixed import circle_intersections

SC = "/tmp/hn"
g = build_G()
field = g.vertices[0].x.field
one = field.rational(1)
xs = [float(v.x) for v in g.vertices]
ys = [float(v.y) for v in g.vertices]
cells = collections.defaultdict(list)
for i, (a, b) in enumerate(zip(xs, ys)):
    cells[(int(a // 2), int(b // 2))].append(i)

t0 = time.time()
have = set(g.vertices)
cand = {}
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
                        cand[q] = cand.get(q, 0) + 1
    if i % 400 == 399:
        print(f"    ... {i+1}/{g.n}, {len(cand)} candidate points  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"{len(cand)} candidate single points  [{time.time()-t0:.0f}s]",
      flush=True)

ranked = sorted(cand, key=lambda q: -cand[q])
best = None
for n, v in enumerate(ranked[:400]):
    w = build_graph(list(g.vertices) + [v])
    piv = w.vertices.index(v)
    deg = len(w.adj[piv])
    if deg < 3:
        continue
    rel = ColourRelations(w, 5)
    if not rel.colourable:
        print(f"  *** G + one point is NOT 5-COLOURABLE -- chi(R^2) >= 6 ***",
              flush=True)
        rel.close()
        break
    pr = pressure(rel, piv)
    c0 = rel.calls
    T, ok = cegar_core(rel, piv, limit=150)
    if ok:
        T = minimise_core(rel, piv, T)
        if best is None or len(T) < best[0]:
            best = (len(T), v, [w.vertices[t] for t in T])
            d2 = collections.Counter(str(v.dist2(w.vertices[t])) for t in T)
            print(f"  point with {deg} neighbours, pressure {pr}: CORE "
                  f"{len(T)} in {rel.calls-c0} solves, d^2 "
                  f"{dict(list(d2.items())[:4])}  [{time.time()-t0:.0f}s]",
                  flush=True)
            pickle.dump((v, best[2]), open(f"{SC}/onepoint_core.pkl", "wb"))
            if len(T) <= 3:
                print(f"  *** CORE OF {len(T)} -- BLOCKABLE ***", flush=True)
                rel.close()
                break
    rel.close()
    if n % 25 == 24:
        print(f"    ... {n+1} points, best {best[0] if best else None}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"best single-point core: {best[0] if best else None}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
