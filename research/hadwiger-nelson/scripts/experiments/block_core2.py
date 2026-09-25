#!/usr/bin/env python3
"""Reproduce the mixed-distance core of two, and try to close it.

The hole-saturated Sa graph -- 397 vertices plus the 332 deep holes it takes
while staying 4-colourable, 729 vertices and 4152 edges -- has a pivot of
degree 19 whose colour is forced onto just **two** targets, sitting on two
different circles.

That is past every ceiling proved in this package. All of them are
same-distance facts: Niven leaves only d^2 = 1/3 among the rational radii,
ramification removes the prime-power circles, and the covering condition caps
capacity at 5. A core of two at mixed distances is subject to none of it,
because its targets are on different circles, so blocking asks only for three
pairwise adjacent images -- an ordinary unit triangle, which the plane has.

Everything here is deterministic, so the pivot is found again rather than
remembered: same base, same holes in the same order, same degree ordering.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_Sa
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.mixed import (conflict_isometries, conflict_rotation_set,
                      count_cross_transversals, deep_holes)
from hn.multispindle import cross_blocks, cross_conflict_graph
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "60"))
MIN_DEG = int(os.environ.get("HN_MINDEG", "4"))
LIMIT = int(os.environ.get("HN_LIMIT", "2500"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad")


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
    t0 = time.time()
    pts = list(build_Sa())
    g = build_graph(pts)
    holes = deep_holes(g, min_degree=MIN_DEG, limit=LIMIT)
    known = set(g.vertices)
    fresh = [(d, x) for d, x in holes if x not in known]
    lo, hi, best_added = 0, len(fresh), []
    while lo <= hi:
        mid = (lo + hi) // 2
        cand = [x for _, x in fresh[:mid]]
        if is_k_colorable(build_graph(pts + cand), K)[0]:
            best_added, lo = cand, mid + 1
        else:
            hi = mid - 1
    g2 = build_graph(pts + best_added)
    print(f"saturated: {g2}  ({len(best_added)} holes)  "
          f"[{time.time() - t0:.0f}s]", flush=True)

    order = sorted(range(g2.n), key=lambda v: -len(g2.adj[v]))[:PIVOTS]
    found = None
    for bp in order:
        core = minimal_core(g2, K, bp)
        if core and len(core) <= 2:
            found = (bp, core)
            break
    if not found:
        print("  no core of two among these pivots", flush=True)
        return
    bp, core = found
    pv = g2.vertices[bp]
    d2s = [pv.dist2(g2.vertices[j]) for j in core]
    print(f"  pivot {bp} (degree {len(g2.adj[bp])}): core {core}", flush=True)
    for j, d in zip(core, d2s):
        print(f"    target {j}: d^2 = {d}", flush=True)
    save_certificate(g2, f"{OUT}/core2_graph.json", k=K,
                     claim=f"mixed-distance forced core of {len(core)} at k={K}")

    field = pv.x.field
    ident = Rotation(field.rational(1), field.zero())
    families = [("rotations", conflict_rotation_set(g2, bp, core)),
                ("with reflections", conflict_isometries(g2, bp, core))]
    for name, fam in families:
        rots = [ident] + list(fam)
        adj = cross_conflict_graph(g2, bp, core, rots)
        deg = max((len(v) for v in adj.values()), default=0)
        esc = count_cross_transversals(g2, bp, core, rots, cap=200000)
        blocked = cross_blocks(g2, bp, core, rots)
        print(f"  {name:17s}: {len(rots):3d} copies, cross-degree {deg:3d}, "
              f"escapes {esc}{'+' if esc >= 200000 else ''}, blocks: {blocked}",
              flush=True)
        if blocked:
            print("  *** BLOCKED -- building the union ***", flush=True)
            about = [r.about(pv) for r in rots]
            union = {q: None for q in g2.vertices}
            for f in about:
                for q in g2.vertices:
                    union[f(q)] = None
            u = build_graph(list(union))
            ok = is_k_colorable(u, K)[0]
            print(f"  union: {u}  {K}-colourable: {ok}", flush=True)
            if not ok:
                save_certificate(u, f"{OUT}/union_no{K}col.json", k=K,
                                 claim=f"not {K}-colourable")
                print(f"  *** certificate saved: chi >= {K + 1} ***",
                      flush=True)
            return


if __name__ == "__main__":
    main()
