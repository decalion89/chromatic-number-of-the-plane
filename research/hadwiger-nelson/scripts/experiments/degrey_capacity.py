#!/usr/bin/env python3
"""What capacity does de Grey's own forced disjunction sit at?

The narrowing runs found forcing at k=4 on his graph and got the core down to
11 targets at one distance. Nothing ever asked the question that decides
whether that core could be blocked at all: *which* distance, and what its
cycle is.

Niven answers it. cos t = 1 - 1/(2 d^2) is rational exactly when d^2 is, and
only 0, +-1/2, +-1 are rational cosines of rational multiples of pi -- so a
rational squared distance closes into a cycle at four values only, and three
of those cycles are even and trap nothing. Only d^2 = 1/3 traps, and it traps
two.

So for each distance where his graph forces, this prints the core size beside
the capacity of that distance. A core over capacity is a disjunction no number
of copies can ever close.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from collections import defaultdict

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import load_certificate
from hn.graph import build_graph
from hn.coloring import is_k_colorable
from hn.spindle import SeparationTest
from hn.transversal import niven_capacity, niven_order

K = int(os.environ.get("HN_K", "4"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "6"))
CORE = os.environ.get(
    "HN_CORE",
    "/tmp/hn/f4_core.json")


def main() -> None:
    # It has to be a graph that *has* k-colourings. de Grey's G is
    # 5-chromatic, so asking whether some 4-colouring separates a pivot from
    # its targets is vacuous -- every query is UNSAT because there are no
    # 4-colourings at all, and the "forcing" it reports means nothing. The
    # 359-vertex core is 4-colourable and genuinely forces.
    pts, _doc = load_certificate(CORE)
    g = build_graph(pts)
    ok = is_k_colorable(g, K)[0]
    print(g, f"k={K}, {K}-colourable: {ok}", flush=True)
    if not ok:
        print("  vacuous: no k-colouring exists, every separation is UNSAT",
              flush=True)
        return
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    for bp in order[:PIVOTS]:
        pv = g.vertices[bp]
        byd = defaultdict(list)
        for j in range(g.n):
            if j == bp:
                continue
            d2 = pv.dist2(g.vertices[j])
            if d2.is_rational():
                byd[d2.c[0]].append(j)
        allt = sorted(j for v in byd.values() for j in v)
        st = SeparationTest(g, K, bp, allt)
        try:
            for d2v in sorted(byd, key=lambda v: -len(byd[v])):
                if len(byd[d2v]) < 2:
                    continue
                t0 = time.time()
                sep, core = st.run(subset=byd[d2v])
                if sep:
                    continue
                cap = niven_capacity(d2v)
                n = niven_order(d2v)
                verdict = ("blockable" if cap and len(core) <= cap
                           else "no finite order: the images form a path"
                           if not cap else f"core {len(core)} over capacity {cap}")
                print(f"  pivot {bp:4d}  d^2={str(d2v):>8}  targets {len(byd[d2v]):4d}"
                      f"  core {len(core):3d}  cycle {n if n else '-'}"
                      f"  capacity {cap}  -> {verdict}  [{time.time()-t0:.0f}s]",
                      flush=True)
        finally:
            st.close()


if __name__ == "__main__":
    main()
