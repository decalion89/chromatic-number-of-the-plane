#!/usr/bin/env python3
"""Hunt for configurations born with a small forced core.

The angular criterion says what a blockable target set looks like.  With the
120-degree rotation at d^2 = 1/3, two copies conflict when they pick the same
target or two targets 120 degrees apart; with three copies and two targets
pigeonhole forces a shared target, hence a conflict, hence a block.  With
three targets the copies can all pick differently and escape.

So the core has to be at most 2, and narrowing a core of 34 was the wrong
strategy -- the floor was structural and guaranteed from the start.

The right one is breadth: look at many configurations and keep whichever is
born with the smallest core, rather than pouring rounds into one that started
wide.  Each configuration is cheap; what was expensive was insisting on the
same one.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G, build_Sa, build_Y
from hn.fast import IntBasis, fast_graph_complete, fast_walk
from hn.generate import hex_ball, unit_vectors
from hn.graph import build_graph
from hn.spindle import (SeparationTest, local_ball, spindle_union_auto,
                        triple_spindle_union)

K = int(os.environ.get("HN_K", "4"))
THIRD = Fraction(1, 3)
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def smallest_core(g, scan=14, maxd2=Fraction(40)):
    """The narrowest forced core over a sample of pivots, at any distance."""
    best = None
    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:scan]:
        p = g.vertices[bp]
        groups = {}
        for j in range(g.n):
            if j == bp:
                continue
            d2 = p.dist2(g.vertices[j])
            if not d2.is_rational():
                continue
            v = d2.c[0]
            if Fraction(1, 4) <= v <= maxd2:
                groups.setdefault(v, []).append(j)
        if not groups:
            continue
        allt = sorted({j for js in groups.values() for j in js})
        st = SeparationTest(g, K, bp, allt)
        try:
            for val in sorted(groups):
                sep, core = st.run(subset=groups[val])
                if sep or not core:
                    continue
                if best is None or len(core) < best[0]:
                    best = (len(core), bp, val, core)
                    if len(core) <= 2:
                        return best
        finally:
            st.close()
    return best


def configurations():
    """A spread of starting shapes, cheapest first."""
    G = build_G()
    centre = min(range(G.n), key=lambda v: G.vertices[v].fx ** 2 + G.vertices[v].fy ** 2)
    yield "de Grey G", G
    yield "de Grey Sa", build_graph(build_Sa())
    yield "de Grey Y", build_graph(build_Y())
    for r in (2.0, 2.5, 3.0, 3.5, 4.0):
        ball, _ = local_ball(G, centre, r)
        if ball.n > 40:
            yield f"G ball r={r}", ball
    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    for steps, rad in ((4, 2.0), (5, 2.5), (6, 2.5), (6, 3.0), (7, 3.0)):
        rows = fast_walk(b, b.rows(U), steps=steps, radius=rad, cap=40000)
        yield f"walk s={steps} R={rad}", fast_graph_complete(b, rows)
    for r in (3, 4, 5):
        yield f"hex ball r={r}", build_graph(hex_ball(r))


def main():
    t0 = time.time()
    print(f"hunting for a forced core <= 2 at k={K}", flush=True)
    best_overall = None
    for name, g in configurations():
        t = time.time()
        try:
            res = smallest_core(g)
        except Exception as e:
            print(f"  {name:<22} {type(e).__name__}: {e}", flush=True)
            continue
        if res is None:
            print(f"  {name:<22} n={g.n:>6}  no forced core  [{time.time()-t:.0f}s]", flush=True)
            continue
        size, pivot, val, core = res
        print(f"  {name:<22} n={g.n:>6}  core={size} at d2={val} (pivot {pivot})  "
              f"[{time.time()-t:.0f}s]", flush=True)
        if best_overall is None or size < best_overall[0]:
            best_overall = (size, name, g, pivot, val, core)
        if size <= 2:
            print(f"\n*** core {size} -- spindleable ***", flush=True)
            spun = (spindle_union_auto(g, pivot, core[0])[0] if size == 1
                    else triple_spindle_union(g, pivot))
            ok, _ = is_k_colorable(spun, K, timeout=7200)
            print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
            if ok is False:
                save_certificate(spun, os.path.join(OUT, f"hunt_k{K}.json"), K,
                                 f"no proper {K}-colouring, from a core of size {size}")
                print("  *** CERTIFICATE SAVED ***", flush=True)
                return
    if best_overall:
        print(f"\nbest over all configurations: core {best_overall[0]} "
              f"in {best_overall[1]} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
