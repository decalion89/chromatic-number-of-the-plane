#!/usr/bin/env python3
"""Pivots that are not vertices.

Every pivot tried in this package was already a vertex of the graph, which is
a restriction the argument never asked for. The pivot is a point whose colour
is being constrained; any point of the plane will do, and the useful ones are
those with the most graph vertices exactly one away, since forcing comes from
how tightly a pivot's own neighbourhood is pinned.

The candidates are exact: the intersections of unit circles about pairs of
vertices, which is where a point can have two neighbours at all. Adding one to
the graph keeps it a unit-distance graph and can only add constraint.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.coloring import is_k_colorable
from hn.degrey import build_S, build_Sa, build_Y
from hn.graph import build_graph
from hn.mixed import deep_holes
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
HOLES = int(os.environ.get("HN_HOLES", "12"))
MIN_DEG = int(os.environ.get("HN_MINDEG", "6"))
LIMIT = int(os.environ.get("HN_LIMIT", "1500"))
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
    pts = list(BUILDERS[BASE]())
    g = build_graph(pts)
    print(f"{BASE}: {g}  k={K}", flush=True)
    t0 = time.time()
    holes = deep_holes(g, min_degree=MIN_DEG, limit=LIMIT)
    known = set(g.vertices)
    fresh = [(d, x) for d, x in holes if x not in known]
    print(f"  {len(holes)} holes, {len(fresh)} of them new points  "
          f"[{time.time() - t0:.0f}s]", flush=True)
    if not fresh:
        return
    print(f"  top new degrees: {[d for d, _ in fresh[:8]]}", flush=True)

    added = [x for _, x in fresh[:HOLES]]
    g2 = build_graph(pts + added)
    ok = is_k_colorable(g2, K)[0]
    print(f"  with {len(added)} holes added: {g2}  {K}-colourable: {ok}",
          flush=True)
    if not ok:
        print("  the augmented graph is already uncolourable -- "
              "forcing would be vacuous", flush=True)
        return
    index = {p: i for i, p in enumerate(g2.vertices)}
    best = None
    for x in added:
        bp = index[x]
        core = minimal_core(g2, K, bp)
        if not core:
            print(f"  hole pivot deg {len(g2.adj[bp]):3d}: separable",
                  flush=True)
            continue
        pv = g2.vertices[bp]
        circles = {pv.dist2(g2.vertices[j]) for j in core}
        print(f"  hole pivot deg {len(g2.adj[bp]):3d}: core {len(core)} over "
              f"{len(circles)} circles  [{time.time() - t0:.0f}s]", flush=True)
        if best is None or len(core) < best:
            best = len(core)
            if best <= 2:
                print("  *** core 2: no capacity ceiling applies ***",
                      flush=True)
                return
    print(f"  best core from a hole pivot: {best}", flush=True)


if __name__ == "__main__":
    main()
