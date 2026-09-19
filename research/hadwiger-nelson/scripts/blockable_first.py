#!/usr/bin/env python3
"""Ask which pairs can be blocked, then whether they are forced.

The search has been running in the expensive order. Finding a minimal forced
core costs many SAT calls per pivot; testing whether it blocks was then cheap
and almost always said no. With two targets the blocking test is 2-SAT and
decides in linear time, so the order can be swapped: enumerate pairs, keep the
ones a large family of copies would close, and only then spend a SAT query
asking whether that pair is forced.

One query, because a SeparationTest carries a selector per target and answers
any subset of them on the same CNF -- so a pair is one call, not a rebuild.

Blocking is also monotone upward for a core of two, so every copy that can be
named is worth including: unlike larger cores, an extra copy never makes
things worse here.
"""
import os
import random
import sys
import time

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
PIVOTS = int(os.environ.get("HN_PIVOTS", "25"))
PAIRS = int(os.environ.get("HN_PAIRS", "4000"))
SEED = int(os.environ.get("HN_SEED", "11"))
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
    rnd = random.Random(SEED)
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:PIVOTS]
    t0 = time.time()
    blockable = forced = 0
    for bp in order:
        cand = [j for j in range(g.n) if j != bp and j not in g.adj[bp]]
        if len(cand) < 2:
            continue
        st = SeparationTest(g, K, bp, cand)
        try:
            for _ in range(PAIRS):
                a, b = rnd.sample(cand, 2)
                fam = conflict_isometries(g, bp, [a, b])
                if not fam:
                    continue
                if not blocks_two_targets(g, bp, [a, b], [ident] + fam):
                    continue
                blockable += 1
                sep, _core = st.run(subset=[a, b])
                if sep:
                    continue
                forced += 1
                print(f"  *** pivot {bp}, pair ({a}, {b}): BLOCKED and "
                      f"FORCED, {len(fam) + 1} copies  "
                      f"[{time.time() - t0:.0f}s] ***", flush=True)
                return
        finally:
            st.close()
        print(f"  pivot {bp:4d} done: {blockable} blockable pairs so far, "
              f"{forced} also forced  [{time.time() - t0:.0f}s]", flush=True)
    print(f"  {blockable} blockable, {forced} of them forced", flush=True)


if __name__ == "__main__":
    main()
