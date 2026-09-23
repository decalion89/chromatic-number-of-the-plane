"""The 60-degree stack, over many centres -- and what the stack does to N(p).

After stacking six rotated copies of G about p, the neighbourhood of p is the
C6-orbit of N_G(p), and the orbit of a single point is a regular hexagon.  So

    N_U(p) = h disjoint hexagons,  h = # distinct angle classes mod 60 in N(p)

and there are NO edges between different hexagons: two points a unit apart are
60 degrees apart, so if one sits at alpha_i + 120k and the other at
alpha_j + 120l then alpha_i - alpha_j = 60 mod 120, which forces the two
hexagons to coincide.  The stack therefore produces the most 2-colourable
neighbourhood there is -- h disjoint even cycles -- with burden 10 * 2^h, and
every geometric obstruction has been designed out of it.

That is a proof about the local structure, not about the graph: the six copies
still have to be simultaneously escape-coloured by ONE colouring, which is much
stronger than anything the local picture sees.  So the stack is worth running,
over many centres, ranked by what it actually buys:

    h small   -- fewer patterns to kill
    |U| large -- the copies genuinely overlap less, so more constraint arrives

A centre where |U| is close to n is one the graph was already symmetric about,
and stacking there buys nothing.
"""
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver, neighbourhood

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
HOWMANY = int(sys.argv[2]) if len(sys.argv) > 2 else 20
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
half = F.rational(Fr(1, 2)); sin60 = F.sqrt(3) * half; cos60 = half
print(f"{NAME} n={n}   [{time.time()-t0:.0f}s]", flush=True)

cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
mids = {}
for i in range(n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-2, -1, 0, 1, 2):
        for dy in (-2, -1, 0, 1, 2):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i: continue
                ex, ey = hx[j] - hx[i], hy[j] - hy[i]
                if abs(ex * ex + ey * ey - 4.0) > 1e-9: continue
                if (g.vertices[j] - g.vertices[i]).norm2() != F.rational(Fr(4)): continue
                k = (round((hx[i] + hx[j]) / 2, 9), round((hy[i] + hy[j]) / 2, 9))
                if k not in mids:
                    mids[k] = Point((g.vertices[i].x + g.vertices[j].x) * half,
                                    (g.vertices[i].y + g.vertices[j].y) * half)
print(f"  {len(mids)} antipodal midpoints, exact in the field   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def rot(p, c, k):
    x, y = p.x - c.x, p.y - c.y
    for _ in range(k % 6):
        x, y = x * cos60 - y * sin60, x * sin60 + y * cos60
    return Point(x + c.x, y + c.y)

def cosets(kx, ky):
    """|N(p)| and the number of distinct angle classes mod 60 degrees"""
    nb = [i for i in range(n) if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
    cls = {round(math.atan2(hy[i]-ky, hx[i]-kx) % (math.pi/3), 7) for i in nb}
    cls = {min(c, round(math.pi/3 - c, 7)) if False else c for c in cls}
    return len(nb), len(cls)

def union_size(cpt):
    seen = set()
    for k in range(6):
        for q in g.vertices:
            r = rot(q, cpt, k)
            seen.add((round(r.fx, 9), round(r.fy, 9)))
    return len(seen)

ranked = []
for (kx, ky), pt in mids.items():
    m, h = cosets(kx, ky)
    if m >= 6: ranked.append((h, m, kx, ky, pt))
ranked.sort(key=lambda r: (r[0], -r[1]))
print(f"  {len(ranked)} midpoints with |N| >= 6; (cosets, |N|) best: "
      f"{[(r[0], r[1]) for r in ranked[:10]]}   [{time.time()-t0:.0f}s]", flush=True)

scored = []
for h, m, kx, ky, pt in ranked[:60]:
    scored.append((h, -union_size(pt), m, kx, ky, pt))
scored.sort()
print(f"  ranked by (cosets, -|union|): "
      f"{[(s[0], -s[1], s[2]) for s in scored[:10]]}   [{time.time()-t0:.0f}s]",
      flush=True)

for ci, (h, negU, m, kx, ky, cpt) in enumerate(scored[:HOWMANY]):
    pts = {}
    for k in range(6):
        for q in g.vertices:
            r = rot(q, cpt, k)
            pts[(round(r.fx, 9), round(r.fy, 9))] = r
    U = build_graph(list(pts.values()))
    nbU = neighbourhood(U, cpt)
    ms = MuSolver(U, 5, budget=None)
    r = ms.at_most_two(nbU) if ms.colourable else None
    tag = "placeable" if r else ("*** BLOCKED AT TWO ***" if r is False else "?")
    print(f"  centre {ci}: h={h} |N|={m}->{len(nbU)} n={U.n} "
          f"col={ms.colourable} at_most_two={r}  {tag}   [{time.time()-t0:.0f}s]",
          flush=True)
    if ms.colourable is False or r is False:
        json.dump({"source": NAME, "copies": 6,
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.coeffs],
                               [[t.numerator, t.denominator] for t in q.y.coeffs]]
                              for q in U.vertices]},
                  open(f"{ROOT}/data/stack6_hit_{ci}.json", "w"))
        print("    written", flush=True)
    ms.close()
