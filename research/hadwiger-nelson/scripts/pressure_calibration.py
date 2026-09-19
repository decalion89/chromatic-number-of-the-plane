"""Does pressure grow with k, or is 2 a ceiling?

A core of size r at p needs pressure(p) >= k - r, so what the method can
reach depends entirely on how high pressure can go. Every graph measured so
far comes back at exactly 2. If that is a property of de Grey's particular
graph, a k=5-level construction should beat it; if it is a ceiling, it is a
new obstruction and it explains everything at once -- at k=4 pressure 2 is
k-2 and permits the core of two his proof uses, while at k=5 the same 2 is
only k-3 and permits nothing smaller than three.

So: measure the maximum over every vertex of every graph in the package, at
every k where the question is not vacuous.
"""
import sys, time, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.certify import load_certificate
from hn.degrey import build_G, build_S, build_Sa, build_Y
from hn.forced import ColourRelations, min_colours_on
from hn.graph import build_graph
from hn.mixed import joint_core_union

R = "/home/user/darwin-50/research/hadwiger-nelson"
cases = [("S", build_S()), ("Sa", build_Sa()), ("Y", build_Y()),
         ("G", build_G()),
         ("joint-core union", build_graph(joint_core_union()))]
pts, _ = load_certificate(f"{R}/certificates/two_orbit_409_no3coloring.json")
cases.append(("two-orbit 409", build_graph(pts)))
pts, _ = load_certificate(f"{R}/certificates/moser_spindle_no3coloring.json")
cases.append(("Moser spindle", build_graph(pts)))

t0 = time.time()
print(f"{'graph':>18} {'n':>6} {'k':>3} {'max pressure':>13} "
      f"{'= k -':>6} {'smallest core':>14}", flush=True)
for name, g in cases:
    g = g if hasattr(g, "vertices") else build_graph(g)
    for k in (3, 4, 5, 6):
        rel = ColourRelations(g, k)
        if not rel.colourable:
            rel.close()
            continue
        hist = collections.Counter()
        for v in range(g.n):
            if g.adj[v]:
                hist[min_colours_on(rel, sorted(g.adj[v]))] += 1
        rel.close()
        mx = max(hist) if hist else 0
        print(f"{name:>18} {g.n:>6} {k:>3} {mx:>13} {k-mx:>6} {k-mx:>14}   "
              f"{dict(sorted(hist.items()))}  [{time.time()-t0:.0f}s]",
              flush=True)
