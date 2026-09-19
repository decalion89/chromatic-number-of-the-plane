"""de Grey's own construction, one level up: the dihedral group applied to G.

Sa = D6 . S about the origin -- 60-degree rotation and reflection in y -- and
it folds 468 images onto 397 points. That folding is its rigidity, and Sa
reaches pressure 3 at four colours, which is k-1, the classical spindle
regime, and it does so with no forced pair anywhere: the forced-different
graph on a pivot's circle is exactly its 30 unit-distance edges, all degrees
2, no odd cycle. So pressure 3 needs ambient rigidity, not an exotic pair.

Applied to G about the origin the same group folds nothing -- 1581*12 = 18966
points exactly -- because G is not centred for it. The centres where it DOES
fold are computable, one per ordered pair of vertices, since a 60-degree
rotation about c takes u to v exactly when c = R v - R^2 u; the best fold 224
points, the same count as the best translations.
"""
import sys, time, collections, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_G, DEGREY_FIELD
from hn.forced import ColourRelations, min_colours_on
from hn.geometry import Point, Rotation
from hn.graph import build_graph

SC = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
g = build_G()
f = DEGREY_FIELD
rot60 = Rotation(f.rational(Fr(1, 2)), f.sqrt(3) * f.rational(Fr(1, 2)))
folds = pickle.load(open(f"{SC}/rotcentres.pkl", "rb"))
best_k, key, (i, j) = folds[0]
# the exact centre: c = R v - R^2 u for the representative pair
Ru = rot60(g.vertices[i])
c = Point(rot60(g.vertices[j]).x - rot60(Ru).x,
          rot60(g.vertices[j]).y - rot60(Ru).y)
about = rot60.about(c)
print(f"G = {g}; centre folds {best_k} points", flush=True)

t0 = time.time()
cur = list(g.vertices)
for step in range(1, 6):
    cur = list(dict.fromkeys(cur + [about(q) for q in cur]))
    w = build_graph(cur)
    rel = ColourRelations(w, 5)
    if not rel.colourable:
        print(f"  *** {step} rotations: {w} NOT 5-COLOURABLE -- "
              f"chi(R^2) >= 6 ***", flush=True)
        rel.close(); break
    near = sorted(range(w.n), key=lambda v: -len(w.adj[v]))[:250]
    hist = collections.Counter(min_colours_on(rel, sorted(w.adj[v]))
                               for v in near)
    rel.close()
    mx = max(hist)
    print(f"  after {step} rotations: {w}  pressures over the 250 "
          f"highest degrees {dict(sorted(hist.items()))}  max {mx}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if mx >= 3:
        print(f"  *** PRESSURE {mx} -- a core of {5-mx} is now possible ***",
              flush=True)
        break
