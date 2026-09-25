"""Three unit circles about a sqrt3-triangle: how many colours can they need?

Rung one is mu_5(N(h)) >= 3, and with N(h) a hexagon that means its two
inscribed sqrt3-triangles cannot BOTH be monochromatic -- so refuting one is
enough.  And there is a purely local sufficient condition for that:

    if T = {a,b,c} is monochromatic in colour gamma, then every point adjacent
    to a, b or c avoids gamma, so A = N(a) u N(b) u N(c) is left with four
    colours.  chi(G[A]) >= 5 is then a contradiction, and T can never be
    monochromatic.

A is the union of three unit circles centred on a triangle of side sqrt3.  A
point of one circle has two neighbours on it (60 degrees) and at most two on
each of the others (two circles meet twice), so the max degree is 6 and Brooks
allows chi up to 6.  With only TWO circles the degree is 4 and Brooks caps chi
at 4, which is useless; three is the first case where five is not excluded.

So generate the configuration -- every point reachable from the triangle's own
centre by unit steps along the three circles -- and compute chi exactly.  Five
would give a local certificate that a sqrt3-triangle is never monochromatic,
which is rung one outright.  Six would be the whole problem.
"""
import sys, time, math
from collections import deque
from pysat.solvers import Solver

t0 = time.time()
R3 = math.sqrt(3.0)
CENTRES = [(1.0, 0.0), (-0.5, R3 / 2), (-0.5, -R3 / 2)]      # pairwise sqrt3
TOL = 1e-9

def on_circles(p):
    return [k for k, c in enumerate(CENTRES)
            if abs(math.hypot(p[0] - c[0], p[1] - c[1]) - 1.0) < 1e-7]

def circle_points_at_unit(p, k):
    """the points of circle k at distance exactly 1 from p"""
    cx, cy = CENTRES[k]
    dx, dy = p[0] - cx, p[1] - cy
    d = math.hypot(dx, dy)
    if d < TOL or d > 2.0 - TOL:
        return []
    a = d / 2.0
    h = math.sqrt(max(0.0, 1.0 - a * a))
    mx, my = cx + dx * a / d, cy + dy * a / d
    ux, uy = -dy / d, dx / d
    return [(mx + h * ux, my + h * uy), (mx - h * ux, my - h * uy)]

def build(cap, seed_angle=None):
    """Seeding at the triangle's centre closes after 13 points: that centre is a
    very special place, lying on all three circles at once, and its orbit under
    the unit steps is finite and bipartite.  A generic seed on one circle lies on
    only that circle and its orbit is the one that matters, since a real graph
    may put points anywhere on the three curves."""
    if seed_angle is None:
        seed = (0.0, 0.0)                   # the triangle's centre is on all three
    else:
        cx, cy = CENTRES[0]
        seed = (cx + math.cos(seed_angle), cy + math.sin(seed_angle))
    seen = {(round(seed[0], 7), round(seed[1], 7)): seed}
    order = [seed]; q = deque([seed])
    while q and len(order) < cap:
        p = q.popleft()
        for k in range(3):
            for r in circle_points_at_unit(p, k):
                key = (round(r[0], 7), round(r[1], 7))
                if key not in seen:
                    seen[key] = r; order.append(r); q.append(r)
    pts = list(seen.values())
    E = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(math.hypot(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]) - 1.0) < 1e-7:
                E.append((i, j))
    return pts, E

def chi(n, E, cap=7):
    for k in range(1, cap + 1):
        X = lambda v, c: 1 + v * k + c
        cnf = [[X(v, c) for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(); s.delete()
        if ok: return k
    return cap + 1

best = 0
for name, ang in [("centre", None), ("0 deg", 0.0), ("17 deg", math.radians(17)),
                  ("40 deg", math.radians(40)), ("73.22 deg", math.radians(73.2218)),
                  ("90 deg", math.radians(90)), ("1 rad", 1.0),
                  ("golden", math.radians(137.507764))]:
    for cap in (60, 200, 600, 1400):
        pts, E = build(cap, ang)
        if not E: continue
        k = chi(len(pts), E)
        deg = [0] * len(pts)
        for a, b in E:
            deg[a] += 1; deg[b] += 1
        print(f"  seed {name:<10s} cap={cap:<5d} n={len(pts):<5d} "
              f"edges={len(E):<6d} maxdeg={max(deg)}  chi = {k}   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        best = max(best, k)
        if k >= 5:
            print("  *** five colours needed on three circles ***", flush=True)
            sys.exit(0)
print(f"\n  best chi over every seed and depth: {best}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
