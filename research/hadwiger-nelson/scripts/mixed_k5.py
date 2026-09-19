#!/usr/bin/env python3
"""Mixed-distance forcing at k = 5, which has never been searched here.

Every k=5 search in this package grouped targets by their distance from the
pivot -- 29930 queries across all 1581 pivots of de Grey's graph, plus balls
to 82357 vertices, all reporting nothing forced. Every one of them was inside
the same-distance regime, which Niven caps at d^2 = 1/3 and capacity 2, and
which twelve call sites enforce by skipping irrational squared distances.

Mixing the distances is what moved k=4: the same graph, the same colour count,
34 same-distance against 3 mixed. It has never been tried at five.

de Grey's G is the right object for it. It is 5-colourable, so the question is
not vacuous the way it is at k=4 where G has no colouring at all, and it is
not 4-colourable, so it is as constrained as a graph can be while still
admitting five colours.
"""
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.coloring import is_k_colorable
from hn.degrey import build_G, build_Sa, build_Y
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "5"))
BASE = os.environ.get("HN_BASE", "G")
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))          # 0 = every pivot
# Targets beyond this squared distance are dropped. Taking every non-neighbour
# asks "is there a colouring where the pivot's colour is unique to it", i.e.
# "is G - p (k-1)-colourable" -- on a k-chromatic graph that is the expensive
# UNSAT direction, and it took seven minutes on the first pivot alone. Distant
# targets are useless for blocking anyway: their circles are enormous and
# their images never meet.
MAXD = float(os.environ.get("HN_MAXD", "12"))
BUILDERS = {"G": build_G, "Sa": build_Sa, "Y": build_Y}


def shrink(st, core, budget=3000):
    best = list(core)
    for order in (list(core), list(reversed(core))):
        cur, spent, changed = list(core), 0, True
        while changed and len(cur) > 1 and spent < budget:
            changed = False
            for t in [x for x in order if x in cur]:
                spent += 1
                if spent >= budget:
                    break
                trial = [x for x in cur if x != t]
                if not st.run(subset=trial)[0]:
                    cur, changed = trial, True
                    break
        if len(cur) < len(best):
            best = cur
    return best


def main() -> None:
    g = BUILDERS[BASE]()
    ok = is_k_colorable(g, K)[0]
    print(f"{BASE}: {g}  k={K}, {K}-colourable: {ok}", flush=True)
    if not ok:
        print("  vacuous: no k-colouring exists, every query is UNSAT",
              flush=True)
        return
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    hist, best, t0 = Counter(), None, time.time()
    for n, bp in enumerate(order):
        nb = set(g.adj[bp]) | {bp}
        pv = g.vertices[bp]
        targets = [j for j in range(g.n) if j not in nb
                   and float(pv.dist2(g.vertices[j])) <= MAXD]
        if len(targets) < 2:
            continue
        st = SeparationTest(g, K, bp, targets)
        try:
            sep, core = st.run(subset=targets)
            if sep:
                hist["separable"] += 1
                continue
            small = shrink(st, core)
        finally:
            st.close()
        hist[len(small)] += 1
        if best is None or len(small) < best[0]:
            best = (len(small), bp, small)
            pv = g.vertices[bp]
            circles = {pv.dist2(g.vertices[j]) for j in small}
            print(f"  pivot {bp:5d}: MIXED CORE {len(small)} over "
                  f"{len(circles)} circles  [{time.time() - t0:.0f}s]",
                  flush=True)
        if n % 100 == 99:
            print(f"    ... {n + 1} pivots, best {best[0] if best else None}, "
                  f"{hist['separable']} separable  [{time.time() - t0:.0f}s]",
                  flush=True)
    print(f"  core sizes: {dict(sorted(hist.items(), key=lambda kv: str(kv[0])))}",
          flush=True)
    print(f"  best at k={K}: {best[0] if best else None}  "
          f"[{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
