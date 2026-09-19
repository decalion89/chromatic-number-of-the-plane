#!/usr/bin/env python3
"""Colour the union of rotated copies directly, instead of arguing about it.

The independent-transversal test is a *sufficient* condition for the union of
copies to be uncolourable: no choice of one target per copy survives, so no
colouring does. It is not necessary. The union carries every edge between
every pair of copies, not just the ones between chosen target images, and it
can fail to colour for reasons the abstraction never sees.

So build the union and hand it to the solver. The rotations are the ones that
create conflicts among the forced targets -- those are the copies with
something to contribute -- and everything is still exact: each copy's points
are field elements, and every edge is confirmed by exact arithmetic.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.coloring import is_k_colorable
from hn.degrey import build_S, build_Sa, build_Y
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.mixed import conflict_rotation_set
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
COPIES = int(os.environ.get("HN_COPIES", "8"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "6"))
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y}


def minimal_core(g, k, bp):
    nb = set(g.adj[bp]) | {bp}
    targets = [j for j in range(g.n) if j not in nb]
    st = SeparationTest(g, k, bp, targets)
    try:
        sep, core = st.run(subset=targets)
        if sep:
            return None
        small = list(core)
        changed = True
        while changed and len(small) > 1:
            changed = False
            for t in list(reversed(small)):
                trial = [x for x in small if x != t]
                if not st.run(subset=trial)[0]:
                    small, changed = trial, True
                    break
        return small
    finally:
        st.close()


def main() -> None:
    g = build_graph(BUILDERS[BASE]())
    field = g.vertices[0].x.field
    print(f"{BASE}: {g}  k={K}", flush=True)
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:PIVOTS]
    for bp in order:
        core = minimal_core(g, K, bp)
        if not core:
            continue
        rots = conflict_rotation_set(g, bp, core)[:COPIES]
        if not rots:
            continue
        p = g.vertices[bp]
        pts = {q: None for q in g.vertices}
        for r in rots:
            about = r.about(p)
            for q in g.vertices:
                pts[about(q)] = None
        t0 = time.time()
        u = build_graph(list(pts))
        t1 = time.time()
        ok = is_k_colorable(u, K)[0]
        print(f"  pivot {bp:4d}: core {len(core)}, {len(rots)} copies -> {u}  "
              f"{K}-colourable: {ok}  [build {t1-t0:.0f}s, sat {time.time()-t1:.0f}s]",
              flush=True)
        if not ok:
            print(f"  *** the union of {len(rots)} copies is not "
                  f"{K}-colourable ***", flush=True)
            return


if __name__ == "__main__":
    main()
