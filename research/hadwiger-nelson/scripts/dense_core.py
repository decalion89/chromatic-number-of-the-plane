#!/usr/bin/env python3
"""Strip the boundary off a dense cyclotomic ball and colour what is left.

A ball grown to depth d is mostly boundary: a point that is a sum of d steps
has its own neighbours at depth d+1, which is not there, so it ends up with
four edges while the few interior points carry close to two hundred. Mean
degree 6.24 over a median of 4 and a maximum of 198.

The k-core is exactly the part that is not boundary. Peeling it leaves the
interior, where every point still has its full complement of unit neighbours,
and that is the only part with any chance of being hard to colour.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.coloring import is_k_colorable
from hn.cyclograph import scaled_graph
from hn.cyclotomic import CycloField, moser_rotation, unit_steps_extended


N = int(os.environ.get("HN_N", "33"))
K = int(os.environ.get("HN_K", "4"))
MMAX = int(os.environ.get("HN_MMAX", "1"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "3"))
RADIUS = float(os.environ.get("HN_RADIUS", "1.2"))
CAP = int(os.environ.get("HN_CAP", "150000"))
CORES = [int(c) for c in os.environ.get("HN_CORES", "8,16,24,32,48").split(",")]


def main() -> None:
    F = CycloField(N)
    steps = unit_steps_extended(F, moser_rotation(F), MMAX)
    t0 = time.time()
    g = scaled_graph(F, steps, ROUNDS, RADIUS, cap=CAP)
    print(f"Q(zeta_{N}), {len(steps)} steps, depth {ROUNDS}, radius {RADIUS}: "
          f"{g}  [{time.time() - t0:.0f}s]", flush=True)
    for c in CORES:
        t0 = time.time()
        core = g.k_core(c)
        if core.n == 0:
            print(f"  {c}-core: empty", flush=True)
            break
        deg = 2 * core.m / core.n
        t1 = time.time()
        ok = is_k_colorable(core, K)[0]
        print(f"  {c}-core: {core}  deg~{deg:6.2f}  {K}-colourable: {ok}"
              f"  [core {t1 - t0:.0f}s, sat {time.time() - t1:.0f}s]", flush=True)
        if not ok:
            print(f"  *** {c}-core is not {K}-colourable ***", flush=True)
            return


if __name__ == "__main__":
    main()
