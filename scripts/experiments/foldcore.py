"""A core of three needs heavy overlap, not two nearly disjoint copies.

The translated unions folded only 222 of 1581 vertices, so the two copies are
essentially independent: colour one, colour the other, and nothing is forced.
For a pivot to carry a core the ambient graph has to constrain its circle, and
that takes a fold.

de Grey's G has an exact 60-degree rotation centre that folds 789 of its 1581
vertices -- half the graph onto itself.  The rotated union is the most
overlapped 5-chromatic non-critical graph available here, and it is where a
core of three, if one exists at five colours, should show first.

The forward core construction gets its full budget rather than twelve steps,
and the search runs over the pivots of the folded part, where the constraint
is.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import cos, sin, pi
from hn.degrey import build_G
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, pressure, cegar_core, is_core, min_colours_on

g = build_G()
fld = g.vertices[0].x.field
half = fld.rational(Fraction(1, 2))
rot = Rotation(half, fld.sqrt(3) * half)

xs = [float(v.x) for v in g.vertices]
ys = [float(v.y) for v in g.vertices]
c60, s60 = cos(pi / 3), sin(pi / 3)


def Rf(x, y):
    return c60 * x - s60 * y, s60 * x + c60 * y


cnt, rep = collections.Counter(), {}
for i in range(g.n):
    ux, uy = Rf(*Rf(xs[i], ys[i]))
    for j in range(0, g.n, 3):
        vx, vy = Rf(xs[j], ys[j])
        key = (round(vx - ux, 7), round(vy - uy, 7))
        cnt[key] += 1
        rep.setdefault(key, (i, j))
top = cnt.most_common(4)
print(f"best rotation centres by multiplicity: {[c for _, c in top]}",
      flush=True)

for key, mult in top:
    i, j = rep[key]
    # c = R v - R^2 u, exactly
    u, v = g.vertices[i], g.vertices[j]
    c = Point(rot(v).x - rot(rot(u)).x, rot(v).y - rot(rot(u)).y)
    turn = rot.about(c)
    img = [turn(p) for p in g.vertices]
    seen, allp = {}, []
    for p in list(g.vertices) + img:
        if p not in seen:
            seen[p] = len(allp)
            allp.append(p)
    fold = 2 * g.n - len(allp)
    if fold < 100:
        continue
    u2 = build_graph(allp)
    deg = u2.degrees()
    print(f"\n  centre from ({i},{j}), float multiplicity {mult}: "
          f"{u2.n} vertices, {u2.m} edges, {fold} folded", flush=True)
    rel = ColourRelations(u2, 5)
    hubs = sorted(range(u2.n), key=lambda x: -deg[x])[:4]
    for p in hubs:
        t = time.time()
        pr = pressure(rel, p)
        T, ok = cegar_core(rel, p, limit=60)
        msg = f"core of {len(T)}" if ok else f"no core in 60 steps"
        print(f"    pivot {p} (deg {deg[p]}): pressure {pr}, {msg}"
              f"  [{time.time()-t:.0f}s]", flush=True)
        if ok:
            print(f"      *** CORE OF {len(T)} AT FIVE COLOURS, verified "
                  f"{is_core(rel, p, T)} ***", flush=True)
    break
