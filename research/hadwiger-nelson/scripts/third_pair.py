#!/usr/bin/env python3
"""Forced pairs with one target on the classical circle, d^2 = 1/3.

The first mixed-distance forced pair did not block, and the reason was this
morning's theorem reappearing. Its two circles have cos t = 0.6044 and
cos t = 1/5 -- rational, but 1/5 is not among Niven's values -- so neither
angle is a rational multiple of pi, both circle graphs are paths, and a 2-SAT
instance dies only through an odd implication chain. Bipartite on both sides
means an escape always exists.

The capacity ceiling genuinely does not apply to mixed pairs. The *need for
oddness* survives the change of regime, and that is a sharper target than
either fact alone: the pair must be mixed, and at least one of its circles
must carry an odd-order angle. Niven leaves exactly one rational radius that
does, d^2 = 1/3 -- the classical circle, where two points are adjacent exactly
when 120 degrees apart.

The mechanism it buys is explicit. Three copies at 120 degrees put a's images
on a triangle, so at most one of them can choose a, and at least two must
choose b. Whether that closes is then a question about b's circle alone.

Cheap to search, because a SeparationTest answers any subset on one CNF: for
each target a at d^2 = 1/3, every pair (a, b) is one query.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G, build_S, build_Sa, build_Y
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.mixed import (blocks_two_targets, conflict_reflections,
                      conflict_rotations)
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))
THIRD = Fraction(1, 3)
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y, "G": build_G}


def load(name):
    if name in BUILDERS:
        g = BUILDERS[name]()
        return g if hasattr(g, "vertices") else build_graph(g)
    pts, _doc = load_certificate(name)
    return build_graph(pts)


def full_pool(g, bp, core, ident):
    pool = {("id",): ident}
    for x in core:
        for y in range(g.n):
            if y == bp:
                continue
            for r in conflict_rotations(g, bp, x, y):
                pool[("rot", r.cos, r.sin)] = r
            for r in conflict_reflections(g, bp, x, y):
                pool[("ref", r.cos, r.sin)] = r
    return list(pool.values())


def main() -> None:
    g = load(BASE)
    field = g.vertices[0].x.field
    ident = Rotation(field.rational(1), field.zero())
    print(f"{BASE}: {g}  k={K}, {K}-colourable: {is_k_colorable(g, K)[0]}",
          flush=True)
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    t0, pivots_with, forced = time.time(), 0, 0
    for n, bp in enumerate(order):
        pv = g.vertices[bp]
        thirds, others = [], []
        for j in range(g.n):
            if j == bp or j in g.adj[bp]:
                continue
            d2 = pv.dist2(g.vertices[j])
            (thirds if d2.is_rational() and d2.c[0] == THIRD
             else others).append(j)
        if not thirds:
            continue
        pivots_with += 1
        cand = thirds + others
        st = SeparationTest(g, K, bp, cand)
        try:
            for a in thirds:
                if not st.run(subset=[a])[0]:
                    continue        # the leg is forced alone: classical, r=1
                for b in others:
                    if not st.run(subset=[b])[0]:
                        continue    # so is this one; the pair is trivial
                    if not st.run(subset=[a, b])[0]:
                        forced += 1
                        pool = full_pool(g, bp, [a, b], ident)
                        blocked = blocks_two_targets(g, bp, [a, b], pool)
                        print(f"  pivot {bp}: FORCED ({a} at 1/3, {b} at "
                              f"{pv.dist2(g.vertices[b])}), {len(pool)} "
                              f"copies, blocks: {blocked}"
                              f"  [{time.time() - t0:.0f}s]", flush=True)
                        if blocked:
                            print("  *** FORCED AND BLOCKED ***", flush=True)
                            return
        finally:
            st.close()
        if n % 20 == 19:
            print(f"    ... {n + 1} pivots, {pivots_with} with a 1/3 target, "
                  f"{forced} forced pairs  [{time.time() - t0:.0f}s]",
                  flush=True)
    print(f"  {pivots_with} pivots had a 1/3 target; {forced} forced pairs",
          flush=True)


if __name__ == "__main__":
    main()
