"""The step de Grey took once, taken again one colour up.

Sa does not force at four colours; Y = Sa u Sb does, where Sb is Sa rotated by
the chord-1/4 rotation.  G does not force at five colours -- measured, on
sqrt3 and on the 60 commonest distances.  So try the same move: G u rho(G).

A rotation of chord c exists over de Grey's field exactly when c(4-c) is a
square in F = Q(sqrt3, sqrt5, sqrt7, sqrt11), since rho = x + iy with
x = (2-c)/2 and y^2 = 1 - x^2 = c(4-c)/4.  Rational c therefore works whenever
c(4-c) is a rational square times one of the squarefree divisors of 1155.

For each such rotation, about the origin, build the union and ask the forcing
question at five colours: a vertex whose colour every 5-colouring repeats on
its sqrt3-sphere.  The cross edges between the copies are what tighten the
colourings, exactly as Sb tightens Sa.
"""
import sys, time, math
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
g = build_G()
base = list(g.vertices)
print(f"G: {len(base)} points, {len(list(g.edges()))} edges  "
      f"[{time.time()-t0:.0f}s]", flush=True)

SQFREE = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]


def rotation_of(c):
    """rho with |1 - rho|^2 = c, if it lives in the field."""
    x = Fr(2 - Fr(c), 2)
    val = Fr(c) * (4 - Fr(c)) / 4          # = y^2
    if val <= 0:
        return None
    num, den = val.numerator, val.denominator
    t = num * den
    for d in SQFREE:
        if t % d:
            continue
        r = math.isqrt(t // d)
        if r * r * d != t:
            continue
        # y = r/den * sqrt(d)
        y = F.rational(Fr(r, den)) * (F.rational(1) if d == 1 else F.sqrt(d))
        try:
            return Rotation(F.rational(x), y)
        except ValueError:
            continue
    return None


cands = []
for num in range(1, 40):
    for den in (1, 2, 3, 4, 6, 8, 9, 12, 16):
        c = Fr(num, den)
        if not (0 < c < 4) or c in [x[0] for x in cands]:
            continue
        r = rotation_of(c)
        if r is not None:
            cands.append((c, r))
print(f"{len(cands)} rotations available over de Grey's field: "
      + ", ".join(str(c) for c, _ in cands[:18]) + " ...", flush=True)


def forcing_hits(pts, cap=5):
    gg = build_graph(pts)
    E = list(gg.edges())
    n = len(pts)
    cells = {}
    for i, q in enumerate(pts):
        cells.setdefault((math.floor(q.fx), math.floor(q.fy)), []).append(i)
    sphere = [[] for _ in range(n)]
    for i, q in enumerate(pts):
        cx, cy = math.floor(q.fx), math.floor(q.fy)
        for a in (-2, -1, 0, 1, 2):
            for b in (-2, -1, 0, 1, 2):
                for j in cells.get((cx + a, cy + b), ()):
                    if j == i:
                        continue
                    if abs((q.fx - pts[j].fx) ** 2
                           + (q.fy - pts[j].fy) ** 2 - 3.0) > 1e-7:
                        continue
                    dx, dy = pts[j].x - pts[i].x, pts[j].y - pts[i].y
                    if dx * dx + dy * dy == 3:
                        sphere[i].append(j)

    def x(v, c):
        return 1 + v * cap + c

    cls = [[x(v, c) for c in range(cap)] for v in range(n)]
    for a, b in E:
        for c in range(cap):
            cls.append([-x(a, c), -x(b, c)])
    out, big = [], 0
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            return "NOT 5-COLOURABLE", len(E), 0
        for u in sorted(range(n), key=lambda v: -len(sphere[v])):
            S = sphere[u]
            if len(S) < 4:
                break
            big = max(big, len(S))
            if not sv.solve(assumptions=[x(u, 0)] + [-x(v, 0) for v in S]):
                out.append((u, len(S)))
                break
    return out, len(E), big


for c, rho in cands:
    pts = list(base)
    seen = {(q.fx.__round__(9), q.fy.__round__(9)) for q in base}
    for q in base:
        r = rho(q)
        k = (round(r.fx, 9), round(r.fy, 9))
        if k not in seen:
            seen.add(k)
            pts.append(r)
    hits, m, big = forcing_hits(pts)
    tag = (hits if isinstance(hits, str)
           else (f"FORCING at vertex {hits[0][0]}, sphere {hits[0][1]}"
                 if hits else f"no forcing (largest sphere {big})"))
    print(f"  chord {str(c):>6}: {len(pts)} points, {m} edges -> {tag}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
