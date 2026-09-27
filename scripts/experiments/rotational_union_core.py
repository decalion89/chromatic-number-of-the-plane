"""The most tightly coupled union available, measured for cores.

The wall wants a 5-chromatic host with the pivot non-critical, substantially
overlapping a pressure gadget. G union rho(G) has all three: it contains G so
it is 5-chromatic; removing any vertex leaves a whole copy so nothing is
critical; and the 60-degree rotation about a centre solved from
c = R v - R^2 u folds 789 of 1581 vertices onto vertices -- half the graph,
against 224 for the best translation. The pressure gadget is already inside G,
since the hexagon offsets read off the witness are exactly 0, theta/2, theta
and their reflections.

Tight coupling is what small cores need: a core has to cover the pivot's free
colours in EVERY colouring, and a loosely coupled graph lets a counterexample
park the free colour anywhere.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G, DEGREY_FIELD
from hn.forced import (ColourRelations, cegar_core, is_core, minimise_core,
                       pressure)
from hn.geometry import Point, Rotation
from hn.graph import build_graph

SC = "/tmp/hn"
g = build_G()
f = DEGREY_FIELD
rot60 = Rotation(f.rational(Fr(1, 2)), f.sqrt(3) * f.rational(Fr(1, 2)))
folds = pickle.load(open(f"{SC}/rotcentres.pkl", "rb"))
best_k, key, (i, j) = folds[0]
Ru = rot60(g.vertices[i])
c = Point(rot60(g.vertices[j]).x - rot60(Ru).x,
          rot60(g.vertices[j]).y - rot60(Ru).y)
about = rot60.about(c)
img = [about(q) for q in g.vertices]
pts = list(dict.fromkeys(list(g.vertices) + img))
folded = 2 * g.n - len(pts)
w = build_graph(pts)
rel = ColourRelations(w, 5)
print(f"W = {w}  5-colourable {rel.colourable}, {folded} folded", flush=True)

t0, best = time.time(), None
gset, iset = set(g.vertices), set(img)
shared = [n for n, v in enumerate(w.vertices) if v in gset and v in iset]
order = sorted(shared, key=lambda v: -len(w.adj[v])) + \
    sorted(set(range(w.n)) - set(shared), key=lambda v: -len(w.adj[v]))
for n, bp in enumerate(order):
    c0 = rel.calls
    T, ok = cegar_core(rel, bp, limit=400)
    if not ok:
        if n % 10 == 9:
            print(f"    ... {n+1} pivots, none closed in 400, best "
                  f"{best[0] if best else None}  [{time.time()-t0:.0f}s]",
                  flush=True)
        continue
    raw = len(T)
    T = minimise_core(rel, bp, T)
    if best is None or len(T) < best[0]:
        best = (len(T), bp, [w.vertices[t] for t in T])
        pv = w.vertices[bp]
        d2 = collections.Counter(str(pv.dist2(w.vertices[t])) for t in T)
        print(f"  pivot {bp} deg {len(w.adj[bp])} pressure {pressure(rel, bp)}"
              f": core {raw} -> {len(T)} in {rel.calls-c0} solves, d^2 "
              f"{dict(list(d2.items())[:4])}  [{time.time()-t0:.0f}s]",
              flush=True)
        pickle.dump((bp, best[2]), open(f"{SC}/rot_core.pkl", "wb"))
        if len(T) <= 3:
            print(f"  *** CORE OF {len(T)} AT FIVE COLOURS -- BLOCKABLE ***",
                  flush=True)
            break
print(f"best rotational-union core: {best[0] if best else None}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
rel.close()
