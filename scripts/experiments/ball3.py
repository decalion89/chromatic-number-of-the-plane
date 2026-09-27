"""The ball has to REFUSE four, or the hexagon 2-colours for free.

A single hexagon at h minimises the burden -- one component, twenty patterns --
and a dense ball around it maximises what constrains them.  Built: 8984 points,
degree 9.76, |N(h)| = 6 with its six edges.  And at_most_two still says
placeable, which is no surprise once the missing question is asked:

    is the ball even 5-CHROMATIC, or does it admit four colours?

If four suffice anywhere in it, every colouring has a spare colour in hand and
the hexagon 2-colours without the rest of the graph noticing.  Density around h
is worthless unless the graph is tight at five, and a ball grown from unit
vectors has no reason to be.

So: test it, and if it is loose, merge in a graph that is not.  five_247_c is
5-chromatic, vertex-critical, and lives in the same field, so a translate of it
overlaps the ball wherever we put it -- giving one graph that both refuses four
AND is dense at the target, which no construction here has had at once.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, deque
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
NAME = "five_247_c.json"
R = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 9000
MERGE = (sys.argv[3] if len(sys.argv) > 3 else "yes") == "yes"
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(P)
vecs = {}
for x, y in g0.edges():
    for a, b in ((x, y), (y, x)):
        v = Point(g0.vertices[b].x - g0.vertices[a].x,
                  g0.vertices[b].y - g0.vertices[a].y)
        vecs[(round(v.fx, 9), round(v.fy, 9))] = v
V = list(vecs.values())
O = Point(F.rational(0), F.rational(0))
pts = {(0.0, 0.0): O}
q = deque([O])
while q and len(pts) < CAP:
    p = q.popleft()
    for v in V:
        r = Point(p.x + v.x, p.y + v.y)
        if r.fx * r.fx + r.fy * r.fy > R * R + 1e-9: continue
        k = (round(r.fx, 9), round(r.fy, 9))
        if k not in pts:
            pts[k] = r; q.append(r)
r3 = F.sqrt(3); half = F.rational(Fr(1, 2))
c60, s60 = half, r3 * half
hexring = None
for v in V:
    ring = []; x, y = v.x, v.y; ok = True
    for _ in range(6):
        k = (round(float(x), 9), round(float(y), 9))
        if k not in pts: ok = False; break
        ring.append(k); x, y = x * c60 - y * s60, x * s60 + y * c60
    if ok and len(set(ring)) == 6: hexring = set(ring); break
if MERGE:
    # a translate of the 5-chromatic graph, placed so its densest vertex sits
    # two units from h: close enough to overlap the ball, far enough that it
    # does not itself touch the unit circle about h
    adj0 = defaultdict(set)
    for x, y in g0.edges(): adj0[x].add(y); adj0[y].add(x)
    v0 = max(range(g0.n), key=lambda v: len(adj0[v]))
    tx = F.rational(2) - g0.vertices[v0].x
    ty = F.rational(0) - g0.vertices[v0].y
    for p in P:
        r = Point(p.x + tx, p.y + ty)
        pts[(round(r.fx, 9), round(r.fy, 9))] = r
keep = {}
for k, p in pts.items():
    dd = p.fx * p.fx + p.fy * p.fy
    if abs(dd - 1.0) < 1e-9 and k not in hexring: continue
    keep[k] = p
G = build_graph(list(keep.values())); n = G.n
E = list(G.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
h = next(i for i in range(n) if (G.vertices[i] - O).norm2() == F.rational(0))
nb = sorted(adj[h])
print(f"ball R={R} merge={MERGE}: n={n} edges={len(E)} degree {2*len(E)/n:.2f}; "
      f"|N(h)|={len(nb)}   [{time.time()-t0:.0f}s]", flush=True)

def colourable(k, budget=None):
    X = lambda v, c: 1 + v * k + c
    cnf = [[X(v, c) for c in range(k)] for v in range(n)]
    for x, y in E:
        for c in range(k):
            cnf.append([-X(x, c), -X(y, c)])
    tri = None
    for u in range(n):
        for v in sorted(adj[u]):
            w = adj[u] & adj[v]
            if w: tri = (u, v, min(w)); break
        if tri: break
    if tri:
        for i, v in enumerate(tri): cnf.append([X(v, i)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if budget is None: r = s.solve()
    else:
        s.conf_budget(budget); r = s.solve_limited()
    s.delete(); return r

r4 = colourable(4, 20_000_000)
print(f"  4-colourable: {r4}   [{time.time()-t0:.0f}s]", flush=True)
r5 = colourable(5, 20_000_000)
print(f"  5-colourable: {r5}   [{time.time()-t0:.0f}s]", flush=True)
if r5 is False:
    print("  *** SIX CHROMATIC ***", flush=True)
    json.dump({"radius": R, "merged": MERGE, "field_generators": list(F.gens),
               "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                           [[t.numerator, t.denominator] for t in q.y.c]]
                          for q in G.vertices]},
              open(f"{ROOT}/data/ball3_hit.json", "w"))
    sys.exit(0)
ms = MuSolver(G, 5, budget=None)
r = ms.at_most_two(nb)
print(f"  at_most_two(hexagon) = {r}   "
      f"({'placeable' if r else '*** RUNG ONE ***'})   [{time.time()-t0:.0f}s]",
      flush=True)
ms.close()
