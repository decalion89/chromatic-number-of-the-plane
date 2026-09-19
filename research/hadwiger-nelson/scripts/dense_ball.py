#!/usr/bin/env python3
"""Push density inside a fixed radius and watch the chromatic number.

The Eisenstein lattice is discrete: a ball of radius 2 holds thirteen points
and that is all there is.  Z[zeta_33] has rank 20 over Z and is dense in the
plane, so the same ball holds as many points as asked for, all exactly
representable and all joined by exact unit steps.

That is a question nothing in this package has asked: not "how far out can a
construction reach" but "how much can be packed into a small disc before it
stops being colourable".  A radius the size of de Grey's graph is not needed
if the density inside a smaller one is unbounded.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.coloring import is_k_colorable
from hn.cyclograph import CycloPoint, build_cyclo_graph, cyclo_generated
from hn.cyclotomic import CycloField, moser_rotation, unit_steps_extended

N = int(os.environ.get("HN_N", "33"))
K = int(os.environ.get("HN_K", "4"))
RADIUS = float(os.environ.get("HN_RADIUS", "2.2"))
CAPS = [int(c) for c in os.environ.get("HN_CAPS", "4000,10000,25000,50000").split(",")]
MMAX = int(os.environ.get("HN_MMAX", "0"))
# Rounds, not a cap, decides the shape. A cap truncates the breadth-first
# frontier, so the set is all of one layer and a sliver of the next, and those
# sliver points connect to almost nothing: 330 generators gave average degree
# 5.77 where 66 gave 7.97. A complete layer inside the radius is what has the
# degree.
ROUNDS = int(os.environ.get("HN_ROUNDS", "6"))


def main() -> None:
    F = CycloField(N)
    rho = moser_rotation(F) if N % 11 == 0 else None
    if MMAX and rho is not None:
        steps = unit_steps_extended(F, rho, MMAX)
    else:
        steps = [tuple(Fraction(x) for x in u) for u in F.unit_steps()]
    gens = [("add", s) for s in steps]
    if not MMAX and rho is not None:
        gens += [("mul", rho), ("mul", F.conj(rho))]
    print(f"Q(zeta_{N}) degree {F.degree}, radius {RADIUS}, k={K}, "
          f"{len(steps)} unit steps", flush=True)
    for cap in CAPS:
        t0 = time.time()
        pts = [p for p in cyclo_generated(F, [CycloPoint(F, F.zero())], gens, ROUNDS,
                                          cap=cap, radius=RADIUS)]
        g = build_cyclo_graph(pts)
        build = time.time() - t0
        t0 = time.time()
        ok = is_k_colorable(g, K)[0]
        print(f"  cap {cap:6d}: {g}  deg~{2 * g.m / max(g.n, 1):5.2f}  "
              f"{K}-colourable: {ok}   [build {build:.0f}s, sat {time.time() - t0:.0f}s]",
              flush=True)
        if not ok:
            print(f"  *** not {K}-colourable inside radius {RADIUS} ***", flush=True)
            return


if __name__ == "__main__":
    main()
