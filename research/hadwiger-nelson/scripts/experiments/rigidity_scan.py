"""How rigid is the tightest union? Measured by forced-same pairs.

A graph is uniquely k-colourable exactly when the forced-same relation has k
classes, so the direct measure of how far a graph is from that ideal is simply
how many pairs it forces to agree. de Grey's G forces none at all -- it is
5-vertex-critical, and a critical vertex sits alone in the kth colour in some
colouring, so nothing through it is forced either way. That is the far end of
the axis from what is wanted.

Unions escape criticality. The best coupling available is the 60-degree
rotation about a centre solved from c = R v - R^2 u, which folds 789 of 1581
vertices onto vertices -- half the graph, against 224 for the best
translation. If ANY pair is forced to agree there, the axis has been moved
along for the first time, and the object has a direction to be pushed in. If
none is, rigidity needs something other than stacking copies.

Harvesting makes the scan cheap: one model with c(u) = 0 and c(v) = 1
certifies that u and v can differ, so a few models clear most candidates and
only survivors need a query each.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G, DEGREY_FIELD
from hn.forced import ColourRelations
from hn.geometry import Point, Rotation
from hn.graph import build_graph

SC = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
g = build_G()
f = DEGREY_FIELD
rot60 = Rotation(f.rational(Fr(1, 2)), f.sqrt(3) * f.rational(Fr(1, 2)))
best_k, key, (i, j) = pickle.load(open(f"{SC}/rotcentres.pkl", "rb"))[0]
Ru = rot60(g.vertices[i])
c = Point(rot60(g.vertices[j]).x - rot60(Ru).x,
          rot60(g.vertices[j]).y - rot60(Ru).y)
about = rot60.about(c)
img = [about(q) for q in g.vertices]
pts = list(dict.fromkeys(list(g.vertices) + img))
w = build_graph(pts)
rel = ColourRelations(w, 5)
print(f"W = {w}  5-colourable {rel.colourable}, {2*g.n-len(pts)} folded",
      flush=True)

random.seed(23)
t0, forced, checked = time.time(), [], 0
gset, iset = set(g.vertices), set(img)
shared = [n for n, v in enumerate(w.vertices) if v in gset and v in iset]
print(f"  {len(shared)} shared vertices -- the tightest part", flush=True)
order = sorted(shared, key=lambda v: -len(w.adj[v]))
for n, u in enumerate(order):
    cand = set(range(w.n)) - w.adj[u] - {u}
    # harvest: a model with c(u)=0, c(v)=1 shows v can differ from u
    for _ in range(6):
        if not cand or not rel._s.solve(assumptions=[rel._x(u, 0)]):
            break
        m = set(rel._s.get_model())
        cand -= {v for v in cand if rel._x(v, 0) not in m}
        rel._s.set_phases([random.choice([1, -1])
                           * rel._x(v, random.randrange(5))
                           for v in random.sample(range(w.n), 100)])
    for v in sorted(cand):
        checked += 1
        if rel.same(u, v):
            forced.append((u, v))
            print(f"  *** FORCED SAME: {u}, {v}  d^2 = "
                  f"{w.vertices[u].dist2(w.vertices[v])}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    if n % 15 == 14:
        print(f"    ... {n+1}/{len(order)} pivots, {checked} survivors "
              f"queried, {len(forced)} forced same  [{time.time()-t0:.0f}s]",
              flush=True)
print(f"forced-same pairs in the tightest union: {len(forced)} over {checked} "
      f"queries  [{time.time()-t0:.0f}s]", flush=True)
rel.close()
