"""Close the loop at four colours: the gadget has pressure 3, so a core of ONE.

Pressure 3 at k = 4 is pressure k-1, which is exactly what a core of size one
needs -- a forced pair, the classical spindle regime. The 127-vertex gadget
measures it. If the pair is there, spindling it is the oldest move in the
subject and the union is a 5-chromatic unit-distance graph built from a single
angle, cos theta = 5/6, rather than transcribed from a paper.

That does not raise the bound -- chi >= 5 is de Grey's -- but it validates the
whole chain end to end on a target whose answer is known: pressure, core,
block, union, certificate.
"""
import sys, time, collections, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.coloring import is_k_colorable
from hn.forced import (ColourRelations, cegar_core, is_core,
                       minimise_core, pressure)
from hn.graph import build_graph
from hn.mixed import three_hexagon_gadget

SC = "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad"
field, pivot, pts = three_hexagon_gadget()
g = build_graph(pts)
piv = g.vertices.index(pivot)
t0 = time.time()
rel = ColourRelations(g, 4)
print(f"{g}  4-colourable {rel.colourable}, pressure {pressure(rel, piv)}",
      flush=True)
best = None
order = [piv] + sorted(range(g.n), key=lambda v: -len(g.adj[v]))
seen = set()
for p in order:
    if p in seen:
        continue
    seen.add(p)
    pr = pressure(rel, p)
    if pr < 3:
        continue
    c0 = rel.calls
    T, ok = cegar_core(rel, p, limit=60)
    if ok:
        T = minimise_core(rel, p, T)
    if ok and (best is None or len(T) < best[0]):
        best = (len(T), p, list(T))
        pv = g.vertices[p]
        d2 = collections.Counter(str(pv.dist2(g.vertices[t])) for t in T)
        print(f"  pivot {p} deg {len(g.adj[p])} pressure {pr}: CORE OF "
              f"{len(T)} in {rel.calls-c0} solves, distances "
              f"{dict(list(d2.items())[:5])}  [{time.time()-t0:.0f}s]",
              flush=True)
        if len(T) <= 3:
            print("  *** a core of one or two: spindle-able ***", flush=True)
            pickle.dump((pts, g.vertices[p], [g.vertices[t] for t in T]),
                        open(f"{SC}/gadget_core.pkl", "wb"))
            break
    if len(seen) % 20 == 0:
        print(f"    ... {len(seen)} pivots, best {best[0] if best else None}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
rel.close()
print(f"best core on the gadget at k=4: {best[0] if best else None}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
