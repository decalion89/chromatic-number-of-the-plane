#!/usr/bin/env python3
"""Mixed-distance disjunctions: the one case no obstruction here covers.

Every bound proved in this package is a *same-distance* bound, and they all
come from the same geometric fact: rotating about the pivot keeps each target
at its own distance, so same-distance targets put every image on one circle,
where a point has at most two neighbours at distance 1. Degree 2 gives a union
of paths and cycles, Niven leaves only d^2 = 1/3 among rational radii,
ramification removes the prime-power circles, and the covering condition caps
capacity at 5.

Measured, and it is the whole story of the search: de Grey's own forced
disjunction at k=4 sits at d^2 = 1/3, cycle C_3, capacity 2 -- with a core of
34. Seventeen times over.

Targets at *different* distances escape all of it. Their images occupy C
circles, a point can meet two on each, the degree bound becomes 2C, and there
is no single cycle left to be capped. Twelve call sites in this package group
targets by distance before testing them, so the mixed case has never been
looked at here at all.

This measures the minimal mixed-distance core, then asks the cross-target
conflict machinery whether any available rotation closes it.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.graph import build_graph
from hn.multispindle import (available_rotation_orders, cross_blocks,
                             cross_conflict_graph, rotation_powers,
                             _rotation_of_order)
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "12"))
CORE = os.environ.get(
    "HN_CORE",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/f4_core.json")


def shrink(st, core, budget=400):
    core = list(core)
    changed = True
    while changed and len(core) > 1 and budget > 0:
        changed = False
        for t in list(core):
            budget -= 1
            if budget <= 0:
                break
            trial = [x for x in core if x != t]
            if not st.run(subset=trial)[0]:
                core, changed = trial, True
                break
    return core


def main() -> None:
    pts, _doc = load_certificate(CORE)
    g = build_graph(pts)
    field = g.vertices[0].x.field
    print(g, f"k={K}, {K}-colourable: {is_k_colorable(g, K)[0]}", flush=True)
    orders = [n for n in available_rotation_orders(field) if n >= 3]
    print(f"  rotation orders available here: {orders}", flush=True)

    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:PIVOTS]:
        nb = set(g.adj[bp]) | {bp}
        targets = [j for j in range(g.n) if j not in nb]
        st = SeparationTest(g, K, bp, targets)
        try:
            t0 = time.time()
            sep, core = st.run(subset=targets)
            if sep:
                continue
            small = shrink(st, core)
        finally:
            st.close()
        pv = g.vertices[bp]
        circles = {pv.dist2(g.vertices[j]) for j in small}
        print(f"  pivot {bp:4d}: mixed core {len(small)} over {len(circles)} "
              f"circles  [{time.time() - t0:.0f}s]", flush=True)
        for n in orders:
            rot = _rotation_of_order(n, field)
            if rot is None:
                continue
            rots = rotation_powers(rot, n)
            adj = cross_conflict_graph(g, bp, small, rots)
            deg = max((len(v) for v in adj.values()), default=0)
            blocked = cross_blocks(g, bp, small, rots)
            print(f"     order {n:3d}: max cross-degree {deg:3d}  "
                  f"blocks: {blocked}", flush=True)
            if blocked:
                print(f"  *** BLOCKED: core {small} closed by order {n} ***",
                      flush=True)
                return


if __name__ == "__main__":
    main()
