#!/usr/bin/env python3
"""The k=3 success, moved up a floor.

What worked at k=3 was not the size of the graph but its tightness. The
triangular lattice is *exactly* 3-chromatic -- the tightest structure there
is at that level -- and adding the unit-triangle centroids gave every pivot a
leg on the classical circle without breaking the colouring. The forced pair
fell out at once, and the six copies closed it.

The analogue one floor up is a base that is exactly 4-chromatic and as tight
as that allows: the lattice spindled onto itself, which is how a 4-chromatic
unit-distance graph is built in the first place. Then the same centroids, and
the same hunt.

Each ingredient is added only while the graph still k-colours, since past
that point every separation query is UNSAT by triviality and the forcing it
reports means nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import chromatic_number, is_k_colorable
from hn.generate import hex_ball, spindled
from hn.geometry import SPINDLE
from hn.graph import build_graph
from hn.mixed import blocks_two_targets, two_orbit_block, unit_triangle_centroids
from hn.multispindle import cross_blocks
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
RADIUS = int(os.environ.get("HN_RADIUS", "3"))
SPINDLES = int(os.environ.get("HN_SPINDLES", "1"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/hn")
THIRD = Fraction(1, 3)


def main() -> None:
    pts = list(hex_ball(RADIUS))
    for i in range(SPINDLES):
        g = build_graph(pts)
        pivot = g.vertices[i % g.n]
        extra = [q for q in spindled(pts, pivot, SPINDLE) if q not in set(pts)]
        cand = pts + extra
        if is_k_colorable(build_graph(cand), K)[0]:
            pts = cand
        else:
            print(f"  spindle {i + 1} would break {K}-colourability, stopping",
                  flush=True)
            break
    g = build_graph(pts)
    print(f"base: {g}  chromatic number {chromatic_number(g)[0]}", flush=True)

    cents = [c for c in unit_triangle_centroids(g) if c not in set(g.vertices)]
    # binary search rather than one test per centroid: a big ball has
    # thousands of them and testing each costs a SAT call, which is what made
    # radius 5 time out while radius 3 finished in seconds.
    lo, hi, keep = 0, len(cents), []
    while lo <= hi:
        mid = (lo + hi) // 2
        if is_k_colorable(build_graph(pts + cents[:mid]), K)[0]:
            keep, lo = cents[:mid], mid + 1
        else:
            hi = mid - 1
    g = build_graph(pts + keep)
    print(f"  + {len(keep)} of {len(cents)} centroids: {g}  "
          f"{K}-colourable: {is_k_colorable(g, K)[0]}", flush=True)

    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    t0, withleg = time.time(), 0
    for n, bp in enumerate(order):
        p = g.vertices[bp]
        legs, others = [], []
        for j in range(g.n):
            if j == bp or j in g.adj[bp]:
                continue
            d2 = p.dist2(g.vertices[j])
            if not d2.is_rational():
                continue
            (legs if d2.c[0] == THIRD else others).append(j)
        if not legs or not others:
            continue
        withleg += 1
        st = SeparationTest(g, K, bp, legs + others)
        try:
            for a in legs:
                if not st.run(subset=[a])[0]:
                    continue        # the leg is forced alone: classical, r=1
                for b in others:
                    if st.run(subset=[a, b])[0]:
                        continue
                    if not st.run(subset=[b])[0]:
                        continue    # so is this one; the pair is trivial
                    copies = two_orbit_block(g, bp, a, b)
                    if copies is None:
                        continue
                    assert blocks_two_targets(g, bp, [a, b], copies)
                    assert cross_blocks(g, bp, [a, b], copies)
                    print(f"  pivot {bp}: FORCED ({a} at 1/3, {b} at "
                          f"{p.dist2(g.vertices[b])}) and BLOCKED"
                          f"  [{time.time() - t0:.0f}s]", flush=True)
                    union = {q: None for q in g.vertices}
                    for r in copies:
                        f = r.about(p)
                        for q in g.vertices:
                            union[f(q)] = None
                    u = build_graph(list(union))
                    good = is_k_colorable(u, K)[0]
                    print(f"  union: {u}  {K}-colourable: {good}"
                          f"  [{time.time() - t0:.0f}s]", flush=True)
                    if not good:
                        save_certificate(u, f"{OUT}/tight_no{K}col.json", k=K,
                                         claim=f"not {K}-colourable by the "
                                         f"two-orbit block")
                        print(f"  *** chi >= {K + 1} ***", flush=True)
                    return
        finally:
            st.close()
        if n % 20 == 19:
            print(f"    ... {n + 1} pivots, {withleg} with a leg"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
    print(f"  {withleg} pivots had a leg, none forced", flush=True)


if __name__ == "__main__":
    main()
