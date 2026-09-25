#!/usr/bin/env python3
"""The most constrained 4-colourable graph reachable from de Grey's.

Forcing at k needs the graph to *have* k-colourings -- on a 5-chromatic graph
every separation query is UNSAT by triviality and means nothing. So the
question lives on 4-colourable graphs, and the cores shrink as those get more
constrained: 34 same-distance on the 359-vertex core, 5 mixed on the same
graph, 3 mixed on Sa.

The tightest such graph is a *maximal* one: 4-colourable, but with nothing
left to add. de Grey's G is 5-chromatic, so peeling the fewest vertices that
make it colour leaves exactly that. What comes back is the largest, densest
graph on which the forcing question is still meaningful.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import find_uncolorable_core, is_k_colorable
from hn.degrey import build_G
from hn.graph import build_graph
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "40"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/maximal4.json")


def peel_to_colourable(g, k):
    """Remove vertices until the graph k-colours, taking one core at a time.

    An uncolourable core is a witness that k colours will not do; deleting one
    of its vertices kills that witness. Repeating until the whole graph
    colours leaves a graph with no room left in it.
    """
    keep = list(range(g.n))
    removed = []
    while True:
        sub = g.induced(keep)
        if is_k_colorable(sub, k)[0]:
            return sub, removed
        core = find_uncolorable_core(sub, k)
        if not core:
            return sub, removed
        # drop the core's highest-degree vertex: it is in the most witnesses
        victim = max(core, key=lambda v: len(sub.adj[v]))
        removed.append(keep[victim])
        keep = [v for i, v in enumerate(keep) if i != victim]
        print(f"    peeled {len(removed)} (core {len(core)}), "
              f"{len(keep)} left", flush=True)


def main() -> None:
    g = build_G()
    print(g, flush=True)
    t0 = time.time()
    sub, removed = peel_to_colourable(g, K)
    print(f"  maximal {K}-colourable: {sub}  after removing {len(removed)}"
          f"  [{time.time() - t0:.0f}s]", flush=True)
    save_certificate(sub, OUT, k=K,
                     claim=f"maximal {K}-colourable subgraph of de Grey's G")
    best = None
    order = sorted(range(sub.n), key=lambda v: -len(sub.adj[v]))[:PIVOTS]
    for bp in order:
        nb = set(sub.adj[bp]) | {bp}
        targets = [j for j in range(sub.n) if j not in nb]
        st = SeparationTest(sub, K, bp, targets)
        try:
            sep, core = st.run(subset=targets)
            if sep:
                continue
            small = list(core)
            changed = True
            while changed and len(small) > 1:
                changed = False
                for t in list(reversed(small)):
                    trial = [x for x in small if x != t]
                    if not st.run(subset=trial)[0]:
                        small, changed = trial, True
                        break
        finally:
            st.close()
        if best is None or len(small) < best[0]:
            best = (len(small), bp)
            pv = sub.vertices[bp]
            circles = {pv.dist2(sub.vertices[j]) for j in small}
            print(f"  pivot {bp:4d}: mixed core {len(small)} over "
                  f"{len(circles)} circles  [{time.time() - t0:.0f}s]",
                  flush=True)
            if len(small) <= 2:
                print("  *** core 2: no capacity ceiling applies ***",
                      flush=True)
                return


if __name__ == "__main__":
    main()
