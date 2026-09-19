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
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.degrey import build_G
from hn.spindle import SeparationTest
from hn.transversal import niven_capacity, niven_order

K = int(os.environ.get("HN_K", "4"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "6"))


def main() -> None:
    g = build_G()
    print(g, f"k={K}", flush=True)
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
