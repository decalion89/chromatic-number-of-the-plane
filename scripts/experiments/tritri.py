"""Seed the three-circle configuration at a triangle, where the odd cycle lives.

A generic orbit of the three circles came back bipartite at every seed and depth
tried, up to 1400 vertices.  That is not the whole configuration: orbits ARE the
connected components, so a generic seed only ever reports on a generic
component, and the odd cycles -- if any -- sit in special ones.

The shortest possible odd cycle in a unit-distance graph is a unit triangle, so
look for one with a vertex on each circle.  Parametrise u by its angle on
circle(a); v is then one of the two points of circle(b) a unit from u, and w one
of the two points a unit from both u and v; requiring w to lie on circle(c) is a
single equation in the single unknown, so solutions are expected and can be
found by bisection on the sign of |w - c| - 1.

Seeding the search there puts an odd cycle in the orbit by construction, and
then the question is how far above 2 the chromatic number of that component
climbs.  Four would make a sqrt3-triangle's monochromatic hypothesis locally
refutable at five colours -- rung one, from geometry alone.
"""
import sys, time, math
from collections import deque
from pysat.solvers import Solver

t0 = time.time()
R3 = math.sqrt(3.0)
CENTRES = [(1.0, 0.0), (-0.5, R3 / 2), (-0.5, -R3 / 2)]

def circle_points_at_unit(p, k):
    cx, cy = CENTRES[k]
    dx, dy = p[0] - cx, p[1] - cy
    d = math.hypot(dx, dy)
    if d < 1e-12 or d > 2.0 - 1e-12: return []
    a = d / 2.0
    h = math.sqrt(max(0.0, 1.0 - a * a))
    mx, my = cx + dx * a / d, cy + dy * a / d
    ux, uy = -dy / d, dx / d
    return [(mx + h * ux, my + h * uy), (mx - h * ux, my - h * uy)]

def third_vertex(u, v, s):
    mx, my = (u[0] + v[0]) / 2, (u[1] + v[1]) / 2
    dx, dy = v[0] - u[0], v[1] - u[1]
    d = math.hypot(dx, dy)
    h = math.sqrt(max(0.0, 1.0 - (d / 2) ** 2))
    return (mx + s * h * (-dy / d), my + s * h * (dx / d))

def residual(alpha, i, j):
    cx, cy = CENTRES[0]
    u = (cx + math.cos(alpha), cy + math.sin(alpha))
    vs = circle_points_at_unit(u, 1)
    if len(vs) < 2: return None, None
    v = vs[i]
    w = third_vertex(u, v, 1 if j == 0 else -1)
    return math.hypot(w[0] - CENTRES[2][0], w[1] - CENTRES[2][1]) - 1.0, (u, v, w)

sols = []
for i in (0, 1):
    for j in (0, 1):
        prev = None
        for n in range(4001):
            a = 2 * math.pi * n / 4000
            r, tri = residual(a, i, j)
            if r is None: prev = None; continue
            if prev is not None and prev[1] * r < 0:
                lo, hi = prev[0], a
                for _ in range(80):
                    mid = (lo + hi) / 2
                    rm, _ = residual(mid, i, j)
                    if rm is None: break
                    if rm * residual(lo, i, j)[0] <= 0: hi = mid
                    else: lo = mid
                rr, tri = residual((lo + hi) / 2, i, j)
                if tri and abs(rr) < 1e-9: sols.append(tri)
            prev = (a, r)
print(f"  unit triangles with one vertex on each circle: {len(sols)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
if not sols:
    print("  none -- the three circles carry no triangle at all", flush=True)
    sys.exit(0)
u, v, w = sols[0]
print(f"    example: {tuple(round(x,6) for x in u)}, "
      f"{tuple(round(x,6) for x in v)}, {tuple(round(x,6) for x in w)}", flush=True)
for a, b in ((u, v), (v, w), (w, u)):
    assert abs(math.hypot(a[0]-b[0], a[1]-b[1]) - 1.0) < 1e-7

def build(seeds, cap):
    seen = {}; order = []; q = deque()
    for s in seeds:
        k = (round(s[0], 7), round(s[1], 7))
        if k not in seen:
            seen[k] = s; order.append(s); q.append(s)
    while q and len(order) < cap:
        p = q.popleft()
        for k in range(3):
            for r in circle_points_at_unit(p, k):
                key = (round(r[0], 7), round(r[1], 7))
                if key not in seen:
                    seen[key] = r; order.append(r); q.append(r)
    pts = list(seen.values()); E = []
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
for cap in (30, 100, 300, 800, 1600, 3000):
    pts, E = build([u, v, w], cap)
    k = chi(len(pts), E)
    deg = [0]*len(pts)
    for a, b in E: deg[a] += 1; deg[b] += 1
    print(f"  triangle-seeded cap={cap:<5d} n={len(pts):<5d} edges={len(E):<6d} "
          f"maxdeg={max(deg)}  chi = {k}   [{time.time()-t0:.0f}s]", flush=True)
    best = max(best, k)
    if k >= 4:
        print("  *** four or more colours on three circles ***", flush=True)
        break
print(f"\n  best: {best}   [{time.time()-t0:.0f}s]", flush=True)
