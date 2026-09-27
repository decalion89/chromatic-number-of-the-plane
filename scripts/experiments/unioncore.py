"""The gap is narrower than it looked: a core of THREE is not excluded.

The pressure theorem says a core of size r at p needs pressure(p) >= k - r.
At five colours a core of three therefore needs pressure >= 2 -- and the
measured pressure is EXACTLY 2 everywhere.  The necessary condition is met.
Blocking a core of three is possible too: 27344 misalignment configurations
do it, and the size is pinned at exactly three from both sides.

So the missing piece is not pressure and not blocking.  It is an actual core
of three at five colours, in a graph where the pivot can be forced at all --
which rules out de Grey's G if that really is vertex-critical, but not a
UNION G u (G + t): remove any vertex and a whole 5-chromatic copy survives, so
no pivot gets a colour of its own and the criticality corollary does not bite.

This builds such unions and runs the forward core construction at k = 5 on
their best pivots.  A core of three would be the piece the whole argument has
been missing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from hn.graph import build_graph, UnitDistanceGraph
from hn.geometry import Point
from hn.forced import ColourRelations, pressure, cegar_core, is_core

pts = degrey.build_G(as_graph=False)
g = build_graph(pts)
print(f"G: {g.n} vertices, {g.m} edges", flush=True)

# translate by the most frequent short difference: |d|^2 = 1/3 folds the most
best, seen = None, {}
for i in range(0, g.n, 7):
    for j in range(g.n):
        if i == j:
            continue
        d = (pts[j].x - pts[i].x, pts[j].y - pts[i].y)
        if d[0] * d[0] + d[1] * d[1] == g.vertices[0].field.rational(
                __import__("fractions").Fraction(1, 3)):
            seen[(str(d[0]), str(d[1]))] = d
print(f"  {len(seen)} distinct shifts at |d|^2 = 1/3", flush=True)

results = []
for key, d in list(seen.items())[:4]:
    shifted = [Point(p.x + d[0], p.y + d[1]) for p in pts]
    allp, idx = [], {}
    for p in pts + shifted:
        if p not in idx:
            idx[p] = len(allp)
            allp.append(p)
    u = build_graph(allp)
    overlap = 2 * g.n - u.n
    deg = u.degrees()
    print(f"\n  union: {u.n} vertices, {u.m} edges, {overlap} shared", flush=True)
    rel = ColourRelations(u, 5)
    hubs = sorted(range(u.n), key=lambda v: -deg[v])[:3]
    for p in hubs:
        t = time.time()
        pr = pressure(rel, p)
        T, ok = cegar_core(rel, p, limit=12)
        print(f"    pivot {p} (degree {deg[p]}): pressure {pr}, "
              f"core {'of ' + str(len(T)) if ok else 'NOT found in 12 steps'}"
              f"  [{time.time()-t:.0f}s]", flush=True)
        if ok and len(T) <= 3:
            print(f"      *** CORE OF {len(T)} AT FIVE COLOURS ***", flush=True)
            print(f"      verified: {is_core(rel, p, T)}", flush=True)
