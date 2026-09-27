"""Stand where the forcing set becomes a core of three.

If S uses every colour in every colouring, then for ANY point p,
T = S minus N(p) is a core of p -- N(p) union T contains S, so it uses every
colour. The core's size is |S| - |S intersect N(p)|, so the pivot does all the
work and the requirement is just

    p adjacent to at least |S| - 3 elements of S.

  |S| = 4: one element. Any of the two exact points one away from it and
           placed generically.
  |S| = 5: two elements. The two circle intersections of any pair less than
           two apart -- these always exist.
  |S| = 6: three elements, which needs them concyclic at radius EXACTLY one,
           and then p is their circumcentre.

One trap, met on the first attempt: p must not itself belong to S. A minimal
forcing set loses its property when any element is dropped, so a circumcentre
that happens to be a forcing vertex takes its own target away and the core
evaporates. Require p outside S.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y
from hn.forced import ColourRelations, is_core, min_colours_on, minimise_core
from hn.geometry import Point
from hn.graph import build_graph
from hn.mixed import circle_intersections

K = int(sys.argv[1]) if len(sys.argv) > 1 else 4


def forcing_set(rel, seed_order, limit=400):
    S = []
    for _ in range(limit):
        rel.calls += 1
        if not rel._s.solve(assumptions=[-rel._x(v, 0) for v in S]):
            return S, True
        m = set(rel._s.get_model())
        cand = [v for v in seed_order if v not in S and rel._x(v, 0) in m]
        if not cand:
            return S, False
        S.append(cand[0])
    return S, False


def shrink(rel, S):
    cur, changed = list(S), True
    while changed and len(cur) > 1:
        changed = False
        for t in list(cur):
            trial = [x for x in cur if x != t]
            if trial and min_colours_on(rel, trial) >= rel.k:
                cur, changed = trial, True
                break
    return cur


def circumcentre(a, b, c):
    two = a.x.field.rational(2)
    d = (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) * two
    if d == 0:
        return None
    aa, bb, cc = (a.x*a.x + a.y*a.y, b.x*b.x + b.y*b.y, c.x*c.x + c.y*c.y)
    return Point((aa*(b.y-c.y) + bb*(c.y-a.y) + cc*(a.y-b.y)) / d,
                 (aa*(c.x-b.x) + bb*(a.x-c.x) + cc*(b.x-a.x)) / d)


g = build_Sa()
g = g if hasattr(g, "vertices") else build_graph(g)
one = g.vertices[0].x.field.rational(1)
t0, best = time.time(), None
import random
random.seed(37)
for attempt in range(4):
    order = list(range(g.n))
    random.shuffle(order)
    rel = ColourRelations(g, K)
    if not rel.colourable:
        rel.close(); break
    S, ok = forcing_set(rel, order)
    if not ok:
        rel.close(); continue
    S = shrink(rel, S)
    rel.close()
    pts = [g.vertices[v] for v in S]
    Sset = set(pts)
    cands = []
    for trio in itertools.combinations(range(len(S)), 3):
        o = circumcentre(*(pts[i] for i in trio))
        if o is not None and o.dist2(pts[trio[0]]) == one and o not in Sset:
            cands.append((o, trio))
    for pair in itertools.combinations(range(len(S)), 2):
        a, b = pts[pair[0]], pts[pair[1]]
        if float(a.dist2(b)) > 3.999:
            continue
        for o in circle_intersections(a, one, b, one):
            if o not in Sset:
                cands.append((o, pair))
    for o, sub in cands:
        w = build_graph(list(g.vertices) + [o])
        piv = w.vertices.index(o)
        nb = set(w.adj[piv])
        Sidx = [w.vertices.index(q) for q in pts]
        T = [x for x in Sidx if x not in nb and x != piv]
        r2 = ColourRelations(w, K)
        core = is_core(r2, piv, T)
        small = minimise_core(r2, piv, T) if core else None
        r2.close()
        if core and (best is None or len(small) < best[0]):
            best = (len(small), len(S), len(sub))
            print(f"  |S| = {len(S)}, pivot adjacent to {len(Sidx) - len(T)} "
                  f"of it, degree {len(nb)}: CORE {len(T)} -> minimised "
                  f"{len(small)}  [{time.time()-t0:.0f}s]", flush=True)
print(f"k={K}: best core via a forcing set: {best[0] if best else None}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
