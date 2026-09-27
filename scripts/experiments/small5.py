"""Measure mu_5 where it has never been measured: the SMALL neighbourhoods.

Every mu scan in this project filtered to |N(p)| >= 8, or >= 10, on the
assumption that carrying five colours takes many points to carry them.  The
grading says the assumption is backwards, and the ceiling at four confirms it.

Blocking p means refuting every proper colouring of N(p) that uses at most
four colours -- each one is a way for p to still be placeable, and the rest of
the graph has to kill it.  N(p) is a union of paths and hexagons, so that count
grows like a constant to the power |N(p)|:

    |N| =  5, one path        a few hundred patterns
    |N| =  6, one hexagon     C6 has 732 proper 4-colourings, times the
                              five choices of which colour goes unused
    |N| = 30, five hexagons   of order 10^14

and the floor is |N| >= 5, because five colours need five points to sit on.
So the cheapest possible blocked point has exactly five neighbours, and they
must be pairwise forced apart -- a rainbow.  That is exactly what the ceiling
at four turned out to be: |N| = 4, the cap attained, four points pairwise
apart, at two of the 803 vertices and nowhere else.

So: mu_5 at every candidate point with 5 <= |N| <= 7, plus the exact pattern
burden of each, which is the number the rest of the graph must refute.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, math
from fractions import Fraction as Fr
from itertools import combinations, product
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver

ROOT = HN_DIR
t0 = time.time()

def burden(nb, nadj, k=5):
    """proper colourings of N(p) in k colours that use at most k-1 of them"""
    m = len(nb); pos = {v: i for i, v in enumerate(nb)}
    pairs = [(pos[u], pos[v]) for u in nb for v in nadj[u] if pos[v] > pos[u]]
    tot = 0
    for c in product(range(k), repeat=m):
        if len(set(c)) >= k: continue
        if any(c[a] == c[b] for a, b in pairs): continue
        tot += 1
    return tot

for name in ("five_247_c.json", "five_tuned_1_1.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    cells = defaultdict(list)
    for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
    cand = {}
    for i in range(n):
        cx, cy = int(hx[i] // 2), int(hy[i] // 2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cells.get((cx + dx, cy + dy), ()):
                    if j <= i: continue
                    ex, ey = hx[j] - hx[i], hy[j] - hy[i]
                    dd = ex * ex + ey * ey
                    if dd <= 1e-12 or dd >= 4.0: continue
                    h = math.sqrt(max(0.0, 1.0 - dd / 4.0)); r = math.sqrt(dd)
                    mx, my = (hx[i] + hx[j]) / 2, (hy[i] + hy[j]) / 2
                    ux, uy = -ey / r, ex / r
                    for s in (+1, -1):
                        cand[(round(mx + s * h * ux, 9), round(my + s * h * uy, 9))] = 1
    targets = []
    for (kx, ky) in cand:
        nb = [i for i in range(n) if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
        if 5 <= len(nb) <= 7: targets.append((kx, ky, nb))
    print(f"\n{name} n={n}: {len(targets)} candidate points with 5 <= |N| <= 7 "
          f"(of {len(cand)})   [{time.time()-t0:.0f}s]", flush=True)
    ms = MuSolver(g, 5, budget=None)
    print(f"  base colouring: {ms.colourable}   [{time.time()-t0:.0f}s]", flush=True)
    got = Counter(); bur = {}
    for kx, ky, nb in targets:
        nadj = {u: [v for v in nb if v != u and
                    abs((hx[u]-hx[v])**2 + (hy[u]-hy[v])**2 - 1.0) < 1e-9] for u in nb}
        v = ms.mu(nb)
        got[(len(nb), v)] += 1
        if len(nb) not in bur: bur[len(nb)] = burden(nb, nadj)
        if v is not None and v >= 3:
            print(f"    mu_5 = {v} at ({kx:.6f},{ky:.6f}) |N|={len(nb)}"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
    ms.close()
    print(f"  mu_5 by (|N|, mu): {dict(sorted(got.items()))}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"  patterns the rest of the graph must refute, by |N|: {bur}", flush=True)
