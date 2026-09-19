#!/usr/bin/env python3
"""The k=5 sweep with the field restriction removed.

The earlier sweep only considered targets whose spindle rotation already lived
in the graph's field.  That filter was a mistake: whether a pair is forced
monochromatic is decided by SAT on the graph alone, and no field enters into
it.  The field matters only afterwards, to write the rotated copy down -- and
`spindle_union_auto` widens it on demand, which is free.

In de Grey's graph the filter was discarding a large share of the candidates.
Sampling sixty pivots, the radicals the skipped distances need are sqrt17 (840
candidates), sqrt13 (516) and sqrt19 (278), against sqrt2's four.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.spindle import SeparationTest, spindle_union_auto, triple_spindle_union

K = int(os.environ.get("HN_K", "5"))
MAXD2 = Fraction(os.environ.get("HN_MAXD2", "60"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def groups_for(g, pivot):
    """All targets grouped by rational squared distance, no field filter."""
    p = g.vertices[pivot]
    out = {}
    for j in range(g.n):
        if j == pivot:
            continue
        d2 = p.dist2(g.vertices[j])
        if not d2.is_rational():
            continue
        v = d2.c[0]
        if v < Fraction(1, 4) or v > MAXD2:
            continue
        out.setdefault(v, []).append(j)
    return out


def main():
    t0 = time.time()
    g = build_G()
    print(f"{g}   K={K}  no field filter", flush=True)
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    total_q = 0
    for n, pivot in enumerate(order, 1):
        grp = groups_for(g, pivot)
        if not grp:
            continue
        allt = sorted({j for js in grp.values() for j in js})
        test = SeparationTest(g, K, pivot, allt)
        try:
            for val in sorted(grp):
                total_q += 1
                sep, core = test.run(subset=grp[val])
                if sep:
                    continue
                print(f"  pivot {pivot} d2={val}: FORCED, core {len(core)}", flush=True)
                try:
                    if len(core) == 1:
                        spun, fld = spindle_union_auto(g, pivot, core[0])
                    elif len(core) == 2 and val == Fraction(1, 3):
                        spun, fld = triple_spindle_union(g, pivot), None
                    else:
                        print(f"    core {len(core)} at d2={val}: no spindle for this shape", flush=True)
                        continue
                except (ValueError, AssertionError) as e:
                    print(f"    spindle failed: {e}", flush=True)
                    continue
                print(f"    spindled: {spun} over {fld}", flush=True)
                s2, _ = is_k_colorable(spun, K, timeout=7200)
                print(f"    {K}-colourable: {s2}", flush=True)
                if s2 is False:
                    print(f"\n*** NON-{K}-COLOURABLE GRAPH ***", flush=True)
                    save_certificate(spun, os.path.join(OUT, f"found_k{K}.json"), K,
                                     f"chi(R^2) >= {K+1}: no proper {K}-colouring")
                    return
        finally:
            test.close()
        if n % 100 == 0:
            print(f"  ... {n} pivots, {total_q} separation queries, none forced "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"nothing forced: {len(order)} pivots, {total_q} queries [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
