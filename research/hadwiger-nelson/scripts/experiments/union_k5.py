#!/usr/bin/env python3
"""Forcing at k=5 needs a 5-chromatic graph that is *not* vertex-critical.

Measured, and it explains a thousand separable pivots: de Grey's G is
5-vertex-critical. Removing any vertex makes it 4-colourable -- checked
directly against the solver, not inferred -- and separating a pivot from every
other vertex at k=5 asks exactly that question. So nothing is ever forced
there, whatever the distances, and the mixed-distance scan reporting 1080
separable pivots was reporting criticality rather than a failure of method.

A union of G with a rotated copy has two critical cores, so removing one
vertex leaves the other. If that union still 5-colours, the question stops
being vacuous and forcing has room to exist. If it does not 5-colour, the
search is over by a different door.

Then the hunt: a forced pair with one leg on d^2 = 1/3, which
`hn.mixed.two_orbit_block` closes with six copies. At k=5 that is chi >= 6.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.geometry import rotation_joining
from hn.graph import build_graph
from hn.mixed import blocks_two_targets, two_orbit_block
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "5"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))
ANGLE = Fraction(os.environ.get("HN_ANGLE", "1/3"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/hn")
THIRD = Fraction(1, 3)


def main() -> None:
    g = build_G()
    field = g.vertices[0].x.field
    pv = g.vertices[0]
    about = rotation_joining(ANGLE, field).about(pv)
    pts = {q: None for q in g.vertices}
    for q in g.vertices:
        pts[about(q)] = None
    t0 = time.time()
    u = build_graph(list(pts))
    print(f"G u rho(G): {u}  [{time.time() - t0:.0f}s]", flush=True)

    t0 = time.time()
    ok = is_k_colorable(u, K)[0]
    print(f"  {K}-colourable: {ok}  [{time.time() - t0:.0f}s]", flush=True)
    if not ok:
        save_certificate(u, f"{OUT}/union_no{K}col.json", k=K,
                         claim=f"not {K}-colourable")
        print(f"  *** the union is not {K}-colourable: chi >= {K + 1} ***",
              flush=True)
        return

    order = sorted(range(u.n), key=lambda v: -len(u.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    t0, withleg, forced = time.time(), 0, 0
    for n, bp in enumerate(order):
        p = u.vertices[bp]
        legs, others = [], []
        for j in range(u.n):
            if j == bp or j in u.adj[bp]:
                continue
            d2 = p.dist2(u.vertices[j])
            if not d2.is_rational():
                continue
            (legs if d2.c[0] == THIRD else others).append(j)
        if not legs or not others:
            continue
        withleg += 1
        st = SeparationTest(u, K, bp, legs + others)
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
                    copies = two_orbit_block(u, bp, a, b)
                    blocked = copies is not None and blocks_two_targets(
                        u, bp, [a, b], copies)
                    print(f"  pivot {bp}: FORCED ({a} at 1/3, {b} at "
                          f"{p.dist2(u.vertices[b])}), blocks: {blocked}"
                          f"  [{time.time() - t0:.0f}s]", flush=True)
                    if blocked:
                        print("  *** FORCED AND BLOCKED ***", flush=True)
                        return
        finally:
            st.close()
        if n % 10 == 9:
            print(f"    ... {n + 1} pivots, {withleg} with a 1/3 leg, "
                  f"{forced} forced  [{time.time() - t0:.0f}s]", flush=True)
    print(f"  {withleg} pivots had a leg; {forced} forced pairs", flush=True)


if __name__ == "__main__":
    main()
