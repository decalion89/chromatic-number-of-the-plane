"""The forcing gadget, isolated: hexagons at a vertex, and their tips.

One hexagon never forces.  Its rhombus tip t = n_i + n_{i+1} is adjacent only
to n_i and n_{i+1}, so at four colours c(t) is already the colour outside
{c(0), c(n_i), c(n_{i+1})} or c(0) itself -- avoiding c(0) costs nothing.  The
forcing in Y must come from the tips CONSTRAINING EACH OTHER, which needs more
than one hexagon.

So build exactly that and nothing else: the origin, m hexagons of unit points
at m different angles, and every rhombus tip they carry.  Then ask the sphere
question as m grows, at four colours to see the effect appear, and at five to
see how far it goes.

The angles come from the rotations of de Grey's field -- a chord c admits one
iff c(4-c) is a square there -- so every configuration is exact.
"""
import sys, time, math
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.geometry import DEGREY_FIELD as F, Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

SQFREE = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]


def rotation_of(c):
    x = Fr(2 - Fr(c), 2)
    val = Fr(c) * (4 - Fr(c)) / 4
    if val <= 0:
        return None
    t = val.numerator * val.denominator
    for d in SQFREE:
        if t % d:
            continue
        r = math.isqrt(t // d)
        if r * r * d != t:
            continue
        y = F.rational(Fr(r, val.denominator)) * (F.rational(1) if d == 1
                                                  else F.sqrt(d))
        try:
            return Rotation(F.rational(x), y)
        except ValueError:
            continue
    return None


rots = []
for num in range(1, 40):
    for den in (1, 2, 3, 4, 6, 8, 9, 12, 16):
        c = Fr(num, den)
        if not (0 < c < 4):
            continue
        if any(c == cc for cc, _ in rots):
            continue
        r = rotation_of(c)
        if r is not None:
            rots.append((c, r))
print(f"{len(rots)} rotations available", flush=True)

half = F.rational(Fr(1, 2))
Z6 = Rotation(half, F.sqrt(3) * half)
ORIGIN = Point(F.zero(), F.zero())
ONE = Point(F.rational(1), F.zero())


def hexagon(rot):
    out, p = [], rot(ONE)
    for _ in range(6):
        out.append(p)
        p = Z6(p)
    return out


def gadget(m):
    pts = [ORIGIN]
    seen = {(0.0, 0.0)}
    for c, rot in rots[:m]:
        hexa = hexagon(rot)
        for k, p in enumerate(hexa):
            q = hexa[(k + 1) % 6]
            tip = Point(p.x + q.x, p.y + q.y)
            for z in (p, tip):
                key = (round(z.fx, 9), round(z.fy, 9))
                if key not in seen:
                    seen.add(key)
                    pts.append(z)
    return pts


def sphere_property(pts, k):
    g = build_graph(pts)
    E = list(g.edges())
    n = len(pts)
    S = [j for j in range(1, n)
         if pts[j].x * pts[j].x + pts[j].y * pts[j].y == 3]

    def x(v, c):
        return 1 + v * k + c

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if not sv.solve():
            return "not even k-colourable", len(E), len(S)
        sep = sv.solve(assumptions=[x(0, 0)] + [-x(v, 0) for v in S])
    return ("FORCES" if not sep else "no forcing"), len(E), len(S)


t0 = time.time()
for m in range(1, len(rots) + 1):
    pts = gadget(m)
    if len(pts) > 2600:
        print(f"  stopping at m = {m}: {len(pts)} points", flush=True)
        break
    r4, e4, s4 = sphere_property(pts, 4)
    r5, e5, s5 = sphere_property(pts, 5)
    print(f"  m = {m:2d}: {len(pts):4d} points, {e4:5d} edges, sphere "
          f"{s4:3d} -> 4 colours: {r4:12s} | 5 colours: {r5}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if r5 == "FORCES":
        print("  *** FIVE-COLOUR FORCING ***", flush=True)
        break
