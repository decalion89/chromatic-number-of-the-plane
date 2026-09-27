"""Where is the knot?  In a saved pre_cdcl instance, take graph balls of growing radius around the
last points inserted and ask whether the induced subgraph is 5-colourable.  A ball that is not
would already be a (smaller) non-5-colourable unit-distance graph -- to be verified exactly."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
from collections import deque
from pysat.solvers import Solver
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
d = json.load(open(sys.argv[1])); LAST = int(sys.argv[2]) if len(sys.argv) > 2 else 20; K = 5
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]
key = lambda x, y: (round(x, 9), round(y, 9))
vk = {key(p.fx, p.fy): i for i, p in enumerate(V)}
n = len(V); adj = [set() for _ in range(n)]
PF = [(p.fx, p.fy) for p in V]; UF = [(u.fx, u.fy) for u in U]       # floats once
for i, (px, py) in enumerate(PF):
    for ux, uy in UF:
        j = vk.get(key(px + ux, py + uy))
        if j is not None and j != i: adj[i].add(j); adj[j].add(i)
print(f"graph: n={n}, m={sum(len(a) for a in adj)//2}", flush=True)
src = list(range(n - LAST, n))
dist = {s: 0 for s in src}; dq = deque(src)
while dq:
    v = dq.popleft()
    if dist[v] >= 8: continue
    for w in adj[v]:
        if w not in dist: dist[w] = dist[v] + 1; dq.append(w)
t0 = time.time()
for R in range(1, 9):
    ball = [v for v, dd in dist.items() if dd <= R]; idx = {v: k for k, v in enumerate(ball)}
    E = [(idx[a], idx[b]) for a in ball for b in adj[a] if b in idx and idx[a] < idx[b]]
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(len(ball))]
    for a, b in E:
        for c in range(K): cl.append([-X(a, c), -X(b, c)])
    s = Solver(name="cadical195", bootstrap_with=cl); s.conf_budget(int(__import__('os').environ.get('KBUD', '5000000')))
    r = s.solve_limited(); st = s.accum_stats(); s.delete()
    print(f"  radius {R}: {len(ball)} vertices, {len(E)} edges: 5-colourable {r}  ({st.get('conflicts', 0)} conflicts)   [{time.time()-t0:.0f}s]", flush=True)
    if r is False:
        json.dump({"ball": ball, "radius": R}, open(sys.argv[1].replace(".json", f"_knot_r{R}.json"), "w")); break
    if len(ball) > 0.8 * n: break
