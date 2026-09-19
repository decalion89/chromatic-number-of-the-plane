#!/usr/bin/env python3
"""Apply the spindle argument one level up: 5 colours on de Grey's graph.

de Grey went from 4 to 5 by spindling a graph that was 4-chromatic.  G itself
is 5-chromatic, so the same move at k=5 -- a pair (or a size-2 disjunction at
d^2 = 1/3) that every 5-colouring of G paints alike -- would give a graph with
no proper 5-colouring, i.e. chi(R^2) >= 6.

Cheap to try, which is the point.  The expensive direction is the k=4 UNSAT;
5-colourings of G are found in milliseconds, so each separation query is fast
and thousands of pivots are affordable.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as FIELD
from hn.geometry import required_radical, rotation_joining
from hn.spindle import SeparationTest, spindle_union, triple_spindle_union

K = int(os.environ.get("HN_K", "5"))
MAXD2 = Fraction(os.environ.get("HN_MAXD2", "40"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
ALLOWED = {1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155}


def spindleable_groups(g, pivot):
    """Targets grouped by squared distance, keeping only distances whose
    spindle rotation exists in this field."""
    p = g.vertices[pivot]
    groups, ok = {}, {}
    for j in range(g.n):
        if j == pivot:
            continue
        d2 = p.dist2(g.vertices[j])
        if not d2.is_rational():
            continue
        val = d2.c[0]
        if val < Fraction(1, 4) or val > MAXD2:
            continue
        good = ok.get(val)
        if good is None:
            good = required_radical(val) in ALLOWED
            if good:
                try:
                    rotation_joining(val, FIELD)
                except ValueError:
                    good = False
            ok[val] = good
        if good:
            groups.setdefault(val, []).append(j)
    return groups


def main():
    t = time.time()
    g = build_G()
    print(f"{g}   [built {time.time()-t:.0f}s]", flush=True)
    sat, _ = is_k_colorable(g, K)
    print(f"G is {K}-colourable: {sat}", flush=True)

    # centre-outwards: the most constrained pivots first
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    checked = 0
    for pivot in order:
        groups = spindleable_groups(g, pivot)
        if not groups:
            continue
        allt = sorted({j for js in groups.values() for j in js})
        test = SeparationTest(g, K, pivot, allt)
        try:
            for val in sorted(groups):
                sep, core = test.run(subset=groups[val])
                if sep:
                    continue
                print(f"  pivot {pivot} d2={val}: FORCED, core size {len(core)}", flush=True)
                usable = len(core) == 1 or (len(core) == 2 and val == Fraction(1, 3))
                if not usable:
                    print(f"    core size {len(core)} at d2={val} not spindleable", flush=True)
                    continue
                spun = (spindle_union(g, pivot, core[0], FIELD) if len(core) == 1
                        else triple_spindle_union(g, pivot))
                print(f"    spindled: {spun}", flush=True)
                s2, _ = is_k_colorable(spun, K, timeout=7200)
                print(f"    {K}-colourable: {s2}", flush=True)
                if s2 is False:
                    print(f"\n*** NON-{K}-COLOURABLE UNIT-DISTANCE GRAPH ***", flush=True)
                    save_certificate(
                        spun, os.path.join(OUT, f"found_k{K}.json"), K,
                        f"chi(R^2) >= {K+1}: no proper {K}-colouring",
                    )
                    return
        finally:
            test.close()
        checked += 1
        if checked % 25 == 0:
            print(f"  ... {checked} pivots checked, none forced [{time.time()-t:.0f}s]", flush=True)
    print(f"no forced structure at k={K} over {checked} pivots [{time.time()-t:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
