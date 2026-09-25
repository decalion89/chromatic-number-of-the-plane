#!/usr/bin/env python3
"""Hunt a mixed-distance forced core of two.

The capacity ceiling is a same-distance fact: targets at one distance put
every image on one circle, degree 2, a cycle, and Niven leaves only d^2 = 1/3
among the rational radii -- capacity 2, against de Grey's own core of 34.

A core of *two* at mixed distances is different in kind. Its two targets sit
on different circles, so blocking asks for three pairwise adjacent images
drawn from two circles, which is an ordinary unit triangle and carries no
ceiling at all. Mixing already takes the cores from 34 to 5.

So this scans pivots for the smallest mixed core it can find, over whatever
4-colourable graphs it is given. Two is the target; three would already be
past everything proved here.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from collections import Counter

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_S, build_Sa, build_Y
from hn.graph import build_graph
from hn.spindle import SeparationTest

# Every one of these is 4-colourable, so the question is not vacuous on them.
# More constraint should mean smaller cores, and Y is twice Sa.
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y}

K = int(os.environ.get("HN_K", "4"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))          # 0 = every pivot
BUDGET = int(os.environ.get("HN_BUDGET", "2000"))
CORE = os.environ.get(
    "HN_CORE",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/f4_core.json")


def shrink(st, core, budget):
    """Greedy, but restarted from every target rather than the first that works.

    A minimal core is not unique and greedy order decides which one is found;
    dropping the *largest-index* target first lands somewhere else than
    dropping the first, and on this problem the difference is the whole game.
    """
    best = list(core)
    for order in (list(core), list(reversed(core))):
        cur = list(core)
        changed, spent = True, 0
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
    if CORE in BUILDERS:
        g = build_graph(BUILDERS[CORE]())
    else:
        pts, _doc = load_certificate(CORE)
        g = build_graph(pts)
    ok = is_k_colorable(g, K)[0]
    print(g, f"k={K}, {K}-colourable: {ok}", flush=True)
    if not ok:
        print("  vacuous: no k-colouring exists", flush=True)
        return
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    hist, best, t0 = Counter(), None, time.time()
    for bp in order:
        nb = set(g.adj[bp]) | {bp}
        targets = [j for j in range(g.n) if j not in nb]
        st = SeparationTest(g, K, bp, targets)
        try:
            sep, core = st.run(subset=targets)
            if sep:
                hist["separable"] += 1
                continue
            small = shrink(st, core, BUDGET)
        finally:
            st.close()
        hist[len(small)] += 1
        if best is None or len(small) < best[0]:
            best = (len(small), bp, small)
            pv = g.vertices[bp]
            circles = {pv.dist2(g.vertices[j]) for j in small}
            print(f"  pivot {bp:4d}: core {len(small)} over {len(circles)} "
                  f"circles  [{time.time() - t0:.0f}s]", flush=True)
            if len(small) <= 2:
                print(f"  *** core {len(small)} -- past the ceiling ***",
                      flush=True)
                return
    print(f"  core sizes: {dict(sorted(hist.items(), key=lambda kv: str(kv[0])))}",
          flush=True)
    print(f"  best: {best[0] if best else None}  [{time.time() - t0:.0f}s]",
          flush=True)


if __name__ == "__main__":
    main()
