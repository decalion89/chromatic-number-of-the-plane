"""mu on the richest 1/sqrt3 rings this project can build: sixteen triangles at once.

The ring is a disjoint union of triangles -- two points of it are a unit apart
exactly at 120 degrees, so a coset of 120 holds at most three points and they
are pairwise adjacent, while different cosets never touch.  So chi of a ring is
exactly 3, and mu = 3 says every triangle on it can be given the SAME three
colours.  mu = 4 says two of them are forced to differ, which no measurement in
this project has ever produced, on any set, at five colours.

Coincidence gets less plausible the more triangles there are, and sigma builds
the biggest rings available: one round on the 803-graph gives 8323 points with
rings of up to fifty, which is up to sixteen triangles that would all have to
agree.  The 803 itself tops out at eighteen-point rings, six triangles, and
reads 3.

Cheap per hub once the base colouring is paid for: mu is a nested chain of
assumption calls on one warm solver, and the chain stops at the first yes.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
r11 = F.sqrt(11); sixth = F.rational(Fr(1, 6))
one = F.rational(1); third = F.rational(Fr(1, 3))
def sig(p, c, conj):
    s = -r11 if conj else r11
    dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + (dx - s * dy) * sixth, c.y + (s * dx + dy) * sixth)
pts = {}
for x, y in d["points"]:
    q = Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
    pts[(round(q.fx, 9), round(q.fy, 9))] = q
G0 = build_graph(list(pts.values()))
for c, u in G0.edges():
    for (a, b) in ((c, u), (u, c)):
        for conj in (False, True):
            p = sig(G0.vertices[b], G0.vertices[a], conj)
            k = (round(p.fx, 9), round(p.fy, 9))
            if k not in pts: pts[k] = p
G = build_graph(list(pts.values())); n = G.n
print(f"sigma round 1 of {NAME}: n={n} edges={sum(1 for _ in G.edges())}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
ring = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if abs(dd - 1/3) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == third:
            ring[i].append(j); ring[j].append(i)
def triangles(h):
    vs = ring[h]; out = []
    for a in range(len(vs)):
        for b in range(a+1, len(vs)):
            if (G.vertices[vs[a]] - G.vertices[vs[b]]).norm2() != one: continue
            for c in range(b+1, len(vs)):
                if (G.vertices[vs[a]] - G.vertices[vs[c]]).norm2() == one and \
                   (G.vertices[vs[b]] - G.vertices[vs[c]]).norm2() == one:
                    out.append((vs[a], vs[b], vs[c]))
    return out
counts = {h: len(triangles(h)) for h in ring if len(ring[h]) >= 12}
best = sorted(counts, key=lambda h: -counts[h])[:120]
print(f"  {len(ring)} hubs; richest rings {sorted((len(ring[h]) for h in best), reverse=True)[:6]}"
      f"; triangle counts {sorted((counts[h] for h in best), reverse=True)[:6]}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
ms = MuSolver(G, 5, budget=None)
print(f"  base colouring {ms.colourable}   [{time.time()-t0:.0f}s]", flush=True)
seen = Counter(); top = None
for h in best:
    v = ms.mu(ring[h])
    seen[v] += 1
    if top is None or (v or 0) > top[0]:
        top = (v, h, len(ring[h]), counts[h])
        print(f"    mu(ring) = {v} at h={h}: ring {len(ring[h])} points, "
              f"{counts[h]} triangles   [{time.time()-t0:.0f}s]", flush=True)
    if v is not None and v >= 4:
        print("    *** THE RING PASSES THREE ***", flush=True)
        json.dump({"source": NAME, "hub": h, "ring": ring[h], "mu": v},
                  open(f"{ROOT}/data/ring_mu4.json", "w"))
        break
print(f"  mu on the ring: {dict(sorted(seen.items(), key=str))}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
