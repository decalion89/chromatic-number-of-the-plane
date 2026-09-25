"""A core of three on the three-hexagon gadget, which is what is actually needed.

Pressure 3 was never the target. The pressure theorem says a core of r needs
pressure >= k - r, so pressure 2 at five colours permits r = 3 and forbids
less; and blocking, counting and misalignment together permit r <= 3 and
nothing more. The two meet at exactly three. So the object to find is a core
of THREE, and pressure 2 is not an obstacle to it -- it is the condition that
makes three the right size.

A core is a pressure measurement in disguise: T is a core of p exactly when
min |c(N(p) union T)| = k, so it can be built forwards by counterexamples,
one solve per point added, instead of shrunk.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools, collections, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.forced import ColourRelations, cegar_core, is_core, pressure
from hn.geometry import Point
from hn.graph import build_graph
from hn.mixed import circle_intersections, three_hexagon_gadget

SC = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
field, pivot, pts = three_hexagon_gadget()
one = field.rational(1)
g0 = build_graph(pts)
piv0 = g0.vertices.index(pivot)
aux = [v for v in range(g0.n) if v != piv0 and v not in g0.adj[piv0]]
have = set(pts)
extra = []
for a, b in itertools.combinations(aux, 2):
    if float(g0.vertices[a].dist2(g0.vertices[b])) > 3.999:
        continue
    for q in circle_intersections(g0.vertices[a], one, g0.vertices[b], one):
        if q not in have:
            have.add(q)
            extra.append(q)

t0 = time.time()
for name, points in (("gadget", pts), ("gadget + level two", pts + extra)):
    g = build_graph(points)
    piv = g.vertices.index(pivot)
    rel = ColourRelations(g, 5)
    print(f"{name}: {g}  5-colourable {rel.colourable}, pressure "
          f"{pressure(rel, piv)}", flush=True)
    if not rel.colourable:
        print("  *** NOT 5-COLOURABLE -- chi(R^2) >= 6 ***", flush=True)
        rel.close()
        break
    c0 = rel.calls
    T, ok = cegar_core(rel, piv, limit=200)
    if ok:
        pv = g.vertices[piv]
        d2 = collections.Counter(str(pv.dist2(g.vertices[t])) for t in T)
        print(f"  CORE OF {len(T)} in {rel.calls-c0} solves, verified "
              f"{is_core(rel, piv, T)}; distances {dict(list(d2.items())[:6])}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((points, pivot, [g.vertices[t] for t in T]),
                    open(f"{SC}/hexcore.pkl", "wb"))
    else:
        print(f"  no core within 200 rounds ({rel.calls-c0} solves)  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    rel.close()
