"""What does a point at the ceiling actually look like?

Everything measured so far reads mu = 2.  Before hunting for mu_5 = 5 it is
worth looking at a point that is KNOWN to sit at its ceiling, and the
vertex-critical 5-chromatic graph hands us 803 of them for free:

    G is 5-chromatic and vertex-critical, so for every vertex p the graph
    H = G - p is 4-colourable, and H + p = G is not.  By definition that is
    mu_4(H, p) = 4 -- exactly.  No solving needed to know the value.

So the 803-graph is a sheet of 803 blocked points, one per deletion, and each
comes with its neighbourhood N(p) laid out in the plane.  Three questions:

  1. How big is a ceiling neighbourhood?  If reaching k colours on a bipartite
     N(p) needs many points, mu_5 = 5 needs a very rich point; if it needs
     few, the difficulty is entirely in the global forcing.
  2. What is its shape?  N(p) is a union of paths and even cycles (60 degrees
     is the only unit-distance angle on a radius-1 circle), so the shape is
     just a multiset of component sizes.
  3. Is the ceiling isolated?  In one fixed H, how many OTHER points of the
     plane also reach 4 -- or does p stand alone?

Question 3 is the one that prices the hunt.  If blocked points come in
families, finding one at five means finding a family; if each is a needle,
the search has no gradient to climb.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver, neighbourhood

ROOT = HN_DIR
t0 = time.time()
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
adj = [set() for _ in range(n)]
for x, y in g.edges():
    adj[x].add(y); adj[y].add(x)
print(f"five_247_c  n={n}  edges={sum(len(a) for a in adj)//2}   [{time.time()-t0:.0f}s]", flush=True)

# ---- 1 and 2: the shape of every ceiling neighbourhood, with no solving ----
def shape(nb):
    """component sizes of the induced subgraph on nb, and whether any cycles"""
    s = set(nb); seen = set(); comps = []; cyc = 0
    for v in nb:
        if v in seen: continue
        stack = [v]; seen.add(v); comp = []; ec = 0
        while stack:
            u = stack.pop(); comp.append(u)
            for w in adj[u] & s:
                ec += 1
                if w not in seen:
                    seen.add(w); stack.append(w)
        if ec // 2 >= len(comp): cyc += 1
        comps.append(len(comp))
    return tuple(sorted(comps, reverse=True)), cyc

degs = Counter(); shapes = Counter(); cycled = 0
for p in range(n):
    nb = sorted(adj[p])
    sh, cy = shape(nb)
    degs[len(nb)] += 1; shapes[sh] += 1; cycled += (cy > 0)
print("  degree of a ceiling point (= |N(p)|, mu_4 = 4 exactly):", flush=True)
for k in sorted(degs): print(f"     |N| = {k:2d}  :  {degs[k]:4d} points", flush=True)
print(f"  neighbourhoods containing a cycle: {cycled} of {n}", flush=True)
print("  commonest component shapes:", flush=True)
for sh, c in shapes.most_common(8): print(f"     {sh}  x{c}", flush=True)
mn = min(degs); print(f"  SMALLEST ceiling neighbourhood: |N| = {mn}", flush=True)

# ---- verify on the smallest one, then ask question 3 there ----
p0 = next(p for p in range(n) if len(adj[p]) == mn)
keep = [i for i in range(n) if i != p0]
H = build_graph([g.vertices[i] for i in keep])
print(f"\n  H = G - p0   n={H.n}   [{time.time()-t0:.0f}s]", flush=True)
ms = MuSolver(H, 4, budget=None)
print(f"  H 4-colourable: {ms.colourable}   [{time.time()-t0:.0f}s]", flush=True)
nb0 = neighbourhood(H, g.vertices[p0])
v = ms.mu(nb0)
print(f"  mu_4(H, p0) = {v}   (|N| = {len(nb0)}, criticality says 4)   [{time.time()-t0:.0f}s]", flush=True)

# ---- 3: how many OTHER points of the plane reach 4 in this same H? ----
# candidate points: every point at distance 1 from at least two vertices, i.e.
# the circle intersections.  Collect them by float key, keep the rich ones.
print(f"\n  scanning the plane around H for other ceiling points   [{time.time()-t0:.0f}s]", flush=True)
cand = defaultdict(int)
hx = [q.fx for q in H.vertices]; hy = [q.fy for q in H.vertices]
cells = defaultdict(list)
for i in range(H.n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
import math
seenpt = {}
for i in range(H.n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i: continue
                ex, ey = hx[j] - hx[i], hy[j] - hy[i]
                d2 = ex * ex + ey * ey
                if d2 <= 1e-12 or d2 >= 4.0: continue
                h = math.sqrt(max(0.0, 1.0 - d2 / 4.0))
                mx, my = (hx[i] + hx[j]) / 2, (hy[i] + hy[j]) / 2
                ux, uy = -ey / math.sqrt(d2), ex / math.sqrt(d2)
                for s in (+1, -1):
                    key = (round(mx + s * h * ux, 9), round(my + s * h * uy, 9))
                    cand[key] += 1
rich = sorted(cand.items(), key=lambda kv: -kv[1])
print(f"  {len(cand)} candidate points; richest float counts: "
      f"{[c for _, c in rich[:12]]}   [{time.time()-t0:.0f}s]", flush=True)
# float count 2*C(m,2) for a point with m exact neighbours -> m from count
hits = Counter()
tested = 0
for (kx, ky), c in rich[:400]:
    # rebuild the point exactly is expensive; use the float neighbourhood
    nb = [i for i in range(H.n) if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
    if len(nb) < 4: continue
    r = ms.mu(nb)
    hits[r] += 1; tested += 1
    if r is not None and r >= 3:
        print(f"    mu_4 = {r} at ({kx:.6f},{ky:.6f})  |N|={len(nb)}   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  tested {tested} candidate points in H: {dict(hits)}   [{time.time()-t0:.0f}s]", flush=True)
