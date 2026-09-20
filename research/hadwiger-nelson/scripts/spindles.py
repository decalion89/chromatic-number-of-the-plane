"""Count the gadget, because the gadget is the unit of constraint.

The triangular lattice is rigid at three colours because every triangle spends
all three and that propagates.  At four colours a triangle leaves one spare,
so something bigger has to do the propagating, and the candidate is the Moser
spindle: seven vertices, 4-chromatic, and the smallest thing that is.

A spindle through a vertex a is a pair of points d, g both at squared distance
3 from a -- the long diagonals of two unit rhombi hinged at a -- with
|d - g| = 1.  That is countable directly, and the prediction is sharp: the
triangular lattice is 3-colourable so it cannot contain one at all, while Sa
should be full of them.

If the count is what drives the correlation, then the unit at five colours is
a 5-critical subgraph -- about five hundred vertices against the spindle's
seven -- and a graph would need to be seventy times larger before it packed
them as densely.  That is a quantitative reason for every negative here, and
it says what scale would be needed rather than that none would do.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph

t0 = time.time()


def lattice(side):
    h = F.sqrt(3) / F.rational(2)
    return [Point(F.rational(i) + F.rational(j) / F.rational(2),
                  h * F.rational(j))
            for i in range(side) for j in range(side)]


def spindles(name, pts):
    n = len(pts)
    zf = [(float(p.x), float(p.y)) for p in pts]
    THREE = F.rational(3)
    ONE = F.rational(1)
    cell = {}
    for i, (a, b) in enumerate(zf):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)

    def near(i, r2, tol=1e-7):
        a, b = zf[i]
        rad = r2 ** 0.5
        out = []
        lo, hi = int((a - rad - 1) // 1), int((a + rad + 1) // 1)
        lo2, hi2 = int((b - rad - 1) // 1), int((b + rad + 1) // 1)
        for cx in range(lo, hi + 1):
            for cy in range(lo2, hi2 + 1):
                for j in cell.get((cx, cy), ()):
                    if j == i:
                        continue
                    if abs((a - zf[j][0]) ** 2 + (b - zf[j][1]) ** 2
                           - r2) < tol:
                        out.append(j)
        return out

    total, hubs = 0, 0
    for a in range(n):
        diag = [d for d in near(a, 3.0)
                if (pts[a].x - pts[d].x) ** 2 + (pts[a].y - pts[d].y) ** 2
                == THREE]
        here = 0
        for x in range(len(diag)):
            for y in range(x + 1, len(diag)):
                d, gg = diag[x], diag[y]
                if abs((zf[d][0] - zf[gg][0]) ** 2
                       + (zf[d][1] - zf[gg][1]) ** 2 - 1) > 1e-7:
                    continue
                if (pts[d].x - pts[gg].x) ** 2 + (pts[d].y - pts[gg].y) ** 2 \
                        == ONE:
                    here += 1
        total += here
        if here:
            hubs += 1
    g = build_graph(pts)
    print(f"  {name}: {n} points, {g.m} edges -> {total} Moser spindles "
          f"through {hubs} hubs, {total/n:.2f} per point  "
          f"[{time.time()-t0:.0f}s]", flush=True)


spindles("triangular lattice (chi = 3)", lattice(20))
spindles("Sa (chi = 4)", build_Sa(F))
spindles("Y (chi = 4)", build_Y(F))
spindles("G (chi = 5)", build_G(F, as_graph=False))
