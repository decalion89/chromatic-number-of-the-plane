"""The triangle-centre closure, and every vertex of it tested for a ring hub.

Adjoining the centre of every unit triangle saturates: from the 803-point graph
it reaches 1851 vertices and 7705 edges in four rounds and then adds nothing,
so the closure is a canonical finite object -- the smallest point set containing
five_247_c and closed under taking the centre of a unit triangle.  It carries
10 312 pairs at 1/sqrt3 and rings of up to thirty points, four times what the
graph it came from had.

Still 5-colourable, and still 5-colourable with every one of those 1/sqrt3 pairs
forbidden as well.  But the hub question is finer than the distance question:
it asks whether ONE vertex cannot avoid its OWN ring, and that was only ever
checked at the sixty richest.  Here it is checked at all of them, which is
affordable because escape is the satisfiable direction.
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
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
pts = {}
for x, y in d["points"]:
    q = Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
    pts[(round(q.fx, 9), round(q.fy, 9))] = q
third = F.rational(Fr(1, 3))
while True:
    G = build_graph(list(pts.values()))
    adj = defaultdict(set)
    for x, y in G.edges():
        adj[x].add(y); adj[y].add(x)
    before = len(pts)
    for u in range(G.n):
        for v in sorted(adj[u]):
            if v <= u: continue
            for w in sorted(adj[u] & adj[v]):
                if w <= v: continue
                c = Point((G.vertices[u].x + G.vertices[v].x + G.vertices[w].x) * third,
                          (G.vertices[u].y + G.vertices[v].y + G.vertices[w].y) * third)
                k = (round(c.fx, 9), round(c.fy, 9))
                if k not in pts: pts[k] = c
    if len(pts) == before: break
G = build_graph(list(pts.values())); n = G.n
E = list(G.edges())
print(f"closure of {NAME}: n={n} edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"source": NAME, "field_generators": list(F.gens),
           "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                       [[t.numerator, t.denominator] for t in q.y.c]]
                      for q in G.vertices]},
          open(f"{ROOT}/data/closure_{NAME}", "w"))
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
ring = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if abs(dd - 1/3) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == third:
            ring[i].append(j); ring[j].append(i)
print(f"  {len(ring)} vertices with a ring; largest {max(len(v) for v in ring.values())}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])
hits = []
order = sorted(ring, key=lambda h: -len(ring[h]))
for t, h in enumerate(order):
    cnf = list(base)
    for v in ring[h]:
        for c in range(K):
            cnf.append([-X(h, c), -X(v, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    good = s.solve(); s.delete()
    if not good:
        print(f"  *** RING HUB h={h}, ring {len(ring[h])} ***   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        hits.append(h); break
    if (t + 1) % 200 == 0:
        print(f"    {t+1}/{len(order)} escape   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  ring hubs: {hits or 'none'}   [{time.time()-t0:.0f}s]", flush=True)
