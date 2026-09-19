#!/usr/bin/env python3
"""Hunt a forced disjunction whose targets sit on an odd magic circle.

`hn.transversal` says a same-distance disjunction is trapped only inside an
odd cycle, that a rational squared distance offers exactly one (the triangle at
d^2 = 1/3, capacity 2) by Niven, and that a multiquadratic field offers exactly
that one too, by n | 24.  Everything this package searched lived inside both
bounds.

The escape is an irrational radius 1/(4 sin^2(pi/n)) with n odd and at least 5,
which is the circumradius of a unit n-gon -- and `unit_polygon` builds that
n-gon in closed form.  So: put the n-gon around a pivot, grow unit-distance
structure over it in Q(zeta_33), and ask whether the pivot's colour is forced
onto the polygon.  A core of at most (n+1)/2 is one the rotations can close.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.cyclograph import CycloPoint, build_cyclo_graph, cyclo_generated
from hn.cyclotomic import CycloField, moser_rotation, unit_polygon
from hn.spindle import SeparationTest
from hn.transversal import trapping_bound

K = int(os.environ.get("HN_K", "4"))
ORDER = int(os.environ.get("HN_ORDER", "11"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "3"))
CAP = int(os.environ.get("HN_CAP", "20000"))


def main() -> None:
    F = CycloField(33)
    poly = unit_polygon(F, ORDER)
    rho = moser_rotation(F)
    units = [("add", tuple(Fraction(x) for x in u)) for u in F.unit_steps()]
    # Rotations about the origin by zeta_ORDER map the polygon to itself, so
    # they grow the structure without moving the targets off their circle.
    gens = units + [("mul", rho), ("mul", F.conj(rho)),
                    ("mul", F.zeta(F.n // ORDER))]
    seeds = [CycloPoint(F, F.zero())] + [CycloPoint(F, v) for v in poly]

    print(f"order {ORDER}: capacity {trapping_bound(ORDER)}, k={K}, "
          f"{len(gens)} generators", flush=True)
    for rounds in range(1, ROUNDS + 1):
        t0 = time.time()
        g = build_cyclo_graph(cyclo_generated(F, seeds, gens, rounds, cap=CAP))
        idx = {p.c: i for i, p in enumerate(g.vertices)}
        pivot = idx[F.zero()]
        targets = [idx[v] for v in poly if v in idx]
        st = SeparationTest(g, K, pivot, targets)
        try:
            sep, core = st.run(subset=targets)
        finally:
            st.close()
        note = "separable" if sep else f"FORCED core={len(core)}"
        print(f"  rounds={rounds}: {g}  {len(targets)} targets  {note}  "
              f"[{time.time() - t0:.0f}s]", flush=True)
        if not sep and core is not None and len(core) <= trapping_bound(ORDER):
            print(f"  *** core {len(core)} is within capacity "
                  f"{trapping_bound(ORDER)} ***", flush=True)
            return


if __name__ == "__main__":
    main()
