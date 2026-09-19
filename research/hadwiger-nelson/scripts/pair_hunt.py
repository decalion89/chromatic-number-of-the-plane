#!/usr/bin/env python3
"""Forced pairs, looked for where they can actually be.

Two facts make the search cheap if they are used in the right order.

If the whole target set separates, no subset of it is forced -- so a pivot
whose minimal core does not exist has no forced pair, and every blockable pair
found there is wasted work. That is what the previous pass was doing: hunting
blockable pairs at pivots whose minimal core was five, where nothing can be
forced.

And a forced pair is a forced set of size two, so it lives inside the forced
structure the core already describes. Testing the C(c,2) pairs of a core of
size c is a handful of queries against the same CNF, not a new search.

So: minimal core per pivot, then every pair inside it, then -- only for the
forced ones -- the 2-SAT blocking test against every isometry the field can
name. Both halves are known to be reachable separately; this looks where they
can coincide.
"""
import os
import sys
import time
from collections import Counter
from itertools import combinations

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_S, build_Sa, build_Y
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.mixed import blocks_two_targets, conflict_isometries
from hn.spindle import SeparationTest

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
PIVOTS = int(os.environ.get("HN_PIVOTS", "0"))
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y}


def load(name):
    if name in BUILDERS:
        return build_graph(BUILDERS[name]())
    pts, _doc = load_certificate(name)
    return build_graph(pts)


def main() -> None:
    g = load(BASE)
    field = g.vertices[0].x.field
    ident = Rotation(field.rational(1), field.zero())
    print(f"{BASE}: {g}  k={K}, {K}-colourable: {is_k_colorable(g, K)[0]}",
          flush=True)
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    if PIVOTS:
        order = order[:PIVOTS]
    hist, pairs, t0 = Counter(), 0, time.time()
    for n, bp in enumerate(order):
        cand = [j for j in range(g.n) if j != bp and j not in g.adj[bp]]
        if len(cand) < 2:
            continue
        st = SeparationTest(g, K, bp, cand)
        try:
            sep, core = st.run(subset=cand)
            if sep:
                hist["separable"] += 1
            else:
                small = list(core)
                changed = True
                while changed and len(small) > 1:
                    changed = False
                    for t in list(reversed(small)):
                        trial = [x for x in small if x != t]
                        if not st.run(subset=trial)[0]:
                            small, changed = trial, True
                            break
                hist[len(small)] += 1
                for a, b in combinations(small, 2):
                    if st.run(subset=[a, b])[0]:
                        continue
                    pairs += 1
                    fam = conflict_isometries(g, bp, [a, b])
                    blocked = bool(fam) and blocks_two_targets(
                        g, bp, [a, b], [ident] + fam)
                    print(f"  pivot {bp:4d}: FORCED PAIR ({a}, {b}), "
                          f"{len(fam)} conflict isometries, blocks: {blocked}"
                          f"  [{time.time() - t0:.0f}s]", flush=True)
                    if blocked:
                        print("  *** FORCED AND BLOCKED ***", flush=True)
                        return
        finally:
            st.close()
        if n % 20 == 19:
            print(f"    ... {n + 1} pivots, {pairs} forced pairs, cores "
                  f"{dict(sorted(hist.items(), key=lambda kv: str(kv[0])))}"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
    print(f"  {pairs} forced pairs over {len(order)} pivots; cores "
          f"{dict(sorted(hist.items(), key=lambda kv: str(kv[0])))}", flush=True)


if __name__ == "__main__":
    main()
