#!/usr/bin/env python3
"""End-to-end validation of the two-orbit block, where forcing is easy.

The construction is proved and the only missing ingredient is a forced pair
with a leg on d^2 = 1/3. At k=4 and k=5 that object has not turned up yet, so
the pipeline has never run to the end -- and a proof that has never produced a
certificate is worth less than one that has.

k=3 is where forcing is cheap. Find a forced pair with a 1/3 leg there, apply
the six copies, and hand the union to the solver: if it comes back
uncolourable, the whole chain is validated on a real object, with a
certificate, by a mechanism nobody has used.

The triangular lattice is 3-chromatic and its unit triangles are everywhere,
so their centroids -- exactly 1/sqrt(3) from three corners each -- hand the
leg over for free.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.generate import hex_ball
from hn.graph import build_graph
from hn.mixed import blocks_two_targets, two_orbit_block, unit_triangle_centroids
from hn.multispindle import cross_blocks
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "3"))
RADIUS = int(os.environ.get("HN_RADIUS", "3"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad")
THIRD = Fraction(1, 3)


def main() -> None:
    pts = list(hex_ball(RADIUS))
    g = build_graph(pts)
    cents = [c for c in unit_triangle_centroids(g) if c not in set(g.vertices)]
    g = build_graph(pts + cents)
    ok = is_k_colorable(g, K)[0]
    print(f"hex ball r={RADIUS} + {len(cents)} centroids: {g}  "
          f"{K}-colourable: {ok}", flush=True)
    if not ok:
        print("  vacuous: no k-colouring, every query is UNSAT", flush=True)
        return

    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    t0, forced = time.time(), 0
    for bp in order:
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
                    forced += 1
                    copies = two_orbit_block(g, bp, a, b)
                    if copies is None:
                        continue
                    if not (blocks_two_targets(g, bp, [a, b], copies)
                            and cross_blocks(g, bp, [a, b], copies)):
                        print(f"  pivot {bp}: forced ({a},{b}) but the six "
                              f"copies do NOT block -- the theorem is wrong",
                              flush=True)
                        return
                    print(f"  pivot {bp}: FORCED ({a} at 1/3, {b} at "
                          f"{p.dist2(g.vertices[b])}) and BLOCKED"
                          f"  [{time.time() - t0:.0f}s]", flush=True)
                    about = [r.about(p) for r in copies]
                    union = {q: None for q in g.vertices}
                    for f in about:
                        for q in g.vertices:
                            union[f(q)] = None
                    u = build_graph(list(union))
                    good = is_k_colorable(u, K)[0]
                    print(f"  union of {len(copies)} copies: {u}  "
                          f"{K}-colourable: {good}  [{time.time() - t0:.0f}s]",
                          flush=True)
                    if not good:
                        save_certificate(u, f"{OUT}/two_orbit_no{K}col.json",
                                         k=K, claim=f"not {K}-colourable, by "
                                         f"the two-orbit block on a forced "
                                         f"pair with a leg at d^2=1/3")
                        print(f"  *** certificate saved: chi >= {K + 1} by "
                              f"the two-orbit construction ***", flush=True)
                    return
        finally:
            st.close()
    print(f"  {forced} forced pairs with a leg, none reached the union",
          flush=True)


if __name__ == "__main__":
    main()
