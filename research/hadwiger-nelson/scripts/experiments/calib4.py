"""Calibrate the core builder where the answer is already known.

At four colours de Grey's Sa demonstrably has small cores -- earlier mixed
searches shrank them to 3 by repeated separation queries. If counterexample
construction finds the same thing quickly, the method is sound and the k=5
difficulty is the graph, not the tool; and the sizes it reports say how
rigid a graph has to be before a core becomes small.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_S, build_Sa, build_Y
from hn.forced import ColourRelations, cegar_core, is_core, pressure
from hn.graph import build_graph

for name, builder in (("Sa", build_Sa), ("Y", build_Y)):
    g = builder()
    g = g if hasattr(g, "vertices") else build_graph(g)
    rel = ColourRelations(g, 4)
    print(f"{name} {g}  4-colourable {rel.colourable}", flush=True)
    if not rel.colourable:
        rel.close(); continue
    t0, sizes, best = time.time(), collections.Counter(), None
    for n, bp in enumerate(sorted(range(g.n), key=lambda v: -len(g.adj[v]))):
        c0 = rel.calls
        T, ok = cegar_core(rel, bp, limit=200)
        if ok:
            sizes[len(T)] += 1
            if best is None or len(T) < best[0]:
                best = (len(T), bp, T)
                pv = g.vertices[bp]
                d2 = [str(pv.dist2(g.vertices[t])) for t in T]
                print(f"  pivot {bp} deg {len(g.adj[bp])} pressure "
                      f"{pressure(rel, bp)}: CORE {len(T)} in {rel.calls-c0} "
                      f"solves at d^2 {d2[:6]}  [{time.time()-t0:.0f}s]",
                      flush=True)
        else:
            sizes["none"] += 1
        if n % 40 == 39:
            print(f"    ... {n+1}/{g.n}, {dict(sorted(sizes.items(), key=lambda kv: str(kv[0])))}"
                  f", best {best[0] if best else None}  [{time.time()-t0:.0f}s]",
                  flush=True)
    print(f"  {name}: sizes {dict(sorted(sizes.items(), key=lambda kv: str(kv[0])))}"
          f"  best {best[0] if best else None}  [{time.time()-t0:.0f}s]",
          flush=True)
    rel.close()
