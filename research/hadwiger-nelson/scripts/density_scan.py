#!/usr/bin/env python3
"""How low an independence ratio do the graphs here reach?

alpha(G)/n bounds the largest density a 1-avoiding set can have, so a single
finite unit-distance graph with alpha/n < 1/5 forces chi_m(R^2) >= 6. The
Moser spindle gives 2/7 = 0.2857; the best published bound on m_1 is about
0.2470, by Fourier methods rather than from a graph.

Unlike the search for a 6-chromatic graph this reports a number for every
graph it is given, so there is always something to improve.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.degrey import build_G, build_S, build_Sa, build_Y
from hn.density import independence_number
from hn.geometry import SPINDLE, eisenstein, origin
from hn.graph import build_graph

SOLVER = os.environ.get("HN_SOLVER", "cd19")


def moser():
    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    return build_graph(rh + [SPINDLE(p) for p in rh])


def main() -> None:
    cases = [("Moser spindle", moser),
             ("S", lambda: build_graph(build_S())),
             ("Sa", lambda: build_graph(build_Sa())),
             ("Y", lambda: build_graph(build_Y())),
             ("de Grey G", lambda: build_graph(build_G()))]
    print(f"{'graph':>14} {'n':>6} {'m':>7} {'alpha':>6} {'ratio':>8} "
          f"{'chi_m >=':>9}  time", flush=True)
    for name, fn in cases:
        g = fn()
        t0 = time.time()
        # start the binary search near the spindle ratio rather than at n/2
        guess = max(1, int(g.n * 0.30))
        a, _ = independence_number(g, lo=max(1, guess - int(g.n * 0.12)),
                                   hi=min(g.n, guess + int(g.n * 0.12)),
                                   solver=SOLVER)
        if a == 0:
            a, _ = independence_number(g, solver=SOLVER)
        print(f"{name:>14} {g.n:>6} {g.m:>7} {a:>6} {a / g.n:>8.4f} "
              f"{g.n / a:>9.4f}  [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
