#!/usr/bin/env python3
"""Search the ring Z[omega, sigma] for a unit-distance graph that is not 4-colourable.

Reproducing chi(R^2) >= 5 from our own generator -- rather than transcribing
published coordinates -- is the gate: the same machinery aimed at 5-colourability
is the attack on the open problem.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, os
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np

from hn.generate import lattice_rotations, unit_vectors_multi
from hn.fast import IntBasis, fast_walk, fast_graph_complete
from hn.coloring import is_k_colorable

ROTS = dict(lattice_rotations(100))


def run(subset, depth, steps, radius, k=4, cap=400000, timeout=600):
    U = unit_vectors_multi([ROTS[d] for d in subset], depth=depth)
    basis = IntBasis.covering(U)
    ur = basis.rows(U)
    t = time.time()
    rows = fast_walk(basis, ur, steps=steps, radius=radius, cap=cap)
    t_gen = time.time() - t
    if len(rows) > cap:
        return dict(status="too_big", n=len(rows), t_gen=t_gen)
    t = time.time()
    try:
        g = fast_graph_complete(basis, rows)
    except OverflowError as e:
        return dict(status="overflow", n=len(rows), msg=str(e))
    core = g.k_core(k)
    t_build = time.time() - t
    if core.n == 0:
        return dict(status="empty_core", n=len(rows), m=g.m, t_gen=t_gen, t_build=t_build)
    t = time.time()
    sat, _ = is_k_colorable(core, k, timeout=timeout)
    return dict(
        status={True: "colorable", False: "UNCOLORABLE", None: "timeout"}[sat],
        n=len(rows), m=g.m, core_n=core.n, core_m=core.m,
        t_gen=round(t_gen, 1), t_build=round(t_build, 1), t_sat=round(time.time() - t, 1),
        D=basis.D,
    )


def main():
    configs = []
    for subset in ([3], [3, 7]):
        for depth in (1, 2):
            for radius in (2.0, 2.5, 3.0, 3.5, 4.0):
                for steps in (4, 5, 6, 7, 8):
                    configs.append((subset, depth, steps, radius))
    for subset in ([3, 7, 19],):
        for radius in (2.0, 2.5, 3.0):
            for steps in (3, 4, 5):
                configs.append((subset, 1, steps, radius))

    print(f"{'rots':>10} {'dep':>3} {'st':>3} {'R':>4} | {'|V|':>8} {'|E|':>9} {'core':>8} | {'result':>12} {'gen':>6} {'bld':>6} {'sat':>7}", flush=True)
    for subset, depth, steps, radius in configs:
        r = run(subset, depth, steps, radius)
        tag = r.get("status")
        print(
            f"{str(subset):>10} {depth:>3} {steps:>3} {radius:>4} | "
            f"{r.get('n','-'):>8} {r.get('m','-'):>9} {r.get('core_n','-'):>8} | "
            f"{tag:>12} {r.get('t_gen','-'):>6} {r.get('t_build','-'):>6} {r.get('t_sat','-'):>7}",
            flush=True,
        )
        if tag == "UNCOLORABLE":
            print("\n*** FOUND A NON-4-COLOURABLE UNIT-DISTANCE GRAPH ***", flush=True)
            print(json.dumps(dict(subset=subset, depth=depth, steps=steps, radius=radius, **r), indent=2), flush=True)
            return


if __name__ == "__main__":
    main()
