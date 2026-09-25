"""Does mu_4 ever reach 4 in anything we can build?  The decisive calibration.

Sa, at 397 points, puts 42 of 600 candidate points at mu_4 = 3 and none at 4.
Gluing raises the count.  The question that prices everything is whether ANY
object this project can construct reaches the ceiling at four colours, because
that is the level where the structure demonstrably exists -- free@4 = 0.00 %,
forcing is abundant, and de Grey closed it.

If mu_4 = 4 turns up at a few thousand points, the mechanism works and we know
what enrichment buys.  If it never turns up even in the densest 4-chromatic
carriers here -- mean degree 18.47, against Sa's 9.94 -- then de Grey's 1581
points are doing something that stacking orbits does not, and that is worth
knowing before any more stacking.

The carriers are the orbit stacks: Sa plus k orbits of glue centres, k = 1, 2
and 10, at 2689, 3463 and 6235 points.  Every one of them is 4-colourable in
under a second, so the cost is the mu chains, not the colouring.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from hn.blocked import MuSolver

t0 = time.time()
F = Field((3, 11, 247))
one = F.rational(Fr(1)); half = F.rational(Fr(1, 2))
rot60 = _rot60(F); r30 = Rotation(F.sqrt(3) * half, half)
Sa = build_Sa(F); g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
BEST = [265, 253, 211, 43, 139, 67, 241, 79, 277, 307]
for k in (1, 2, 10):
    cs = []
    for b in BEST[:k]:
        for w in orbit(Sa[b]):
            if w not in cs: cs.append(w)
    seen, U = set(Sa), list(Sa)
    for w in cs:
        rot = g60.about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {k} orbits: n={n} m={m} deg={2.0*m/n:.2f}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    S = set(g.vertices)
    deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
    cand = set()
    for c in deg_order[:40]:
        for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
            for p in g.vertices:
                z = rot(p)
                if z not in S:
                    cand.add(z)
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
    gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
    scored, seenn = [], set()
    for z in cand:
        zx, zy = z.fx, z.fy
        cx, cy = int(zx // 1), int(zy // 1)
        nb = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in cells.get((cx + dx, cy + dy), ()):
                    ex, ey = gx[i] - zx, gy[i] - zy
                    if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                       (g.vertices[i] - z).norm2() == one:
                        nb.append(i)
        if len(nb) >= 8:
            t = tuple(sorted(nb))
            if t not in seenn:
                seenn.add(t); scored.append((len(nb), t))
    scored.sort(key=lambda u: -u[0])
    print(f"    {len(scored)} neighbourhoods of size >= 8, largest "
          f"{scored[0][0] if scored else 0}   [{time.time()-t0:.0f}s]",
          flush=True)
    if not scored:
        continue
    hist = Counter()
    with MuSolver(g, k=4, budget=5_000_000) as ms:
        if not ms.colourable:
            print("    (refuses four already)", flush=True); continue
        for sz, nb in scored[:150]:
            v = ms.mu(list(nb))
            hist[v] += 1
            if v is not None and v >= 4:
                print(f"    *** mu_4 = 4 with |N|={sz}: this carrier + p is "
                      f"5-chromatic on {n+1} vertices ***", flush=True)
    print(f"    mu_4 over the richest {min(150, len(scored))}: {dict(hist)}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
