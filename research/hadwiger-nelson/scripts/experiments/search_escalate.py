#!/usr/bin/env python3
"""Escalate before searching: tighten the graph, then look for forcing at k=5.

de Grey did not hunt for forced pairs inside the Moser spindle.  He built a
far larger and more constrained object -- 20425 vertices from Minkowski sums
and rotated copies -- and found the forcing there.  Our G is 5-chromatic only
barely, so its 5-colourings have room to move and nothing is forced.

This builds unions G u rho(G) about a central pivot, which stay unit-distance
graphs and are markedly more constrained, then runs the separation test on
them.  A 5-colouring of a tighter graph is still found in milliseconds, so the
queries stay cheap even as the vertex count grows.
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
from hn.degrey import build_G
from hn.spindle import SeparationTest, spindle_union_auto

K = int(os.environ.get("HN_K", "5"))
MAXD2 = Fraction(os.environ.get("HN_MAXD2", "40"))
NPIVOTS = int(os.environ.get("HN_PIVOTS", "120"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def groups_for(g, pivot):
    p = g.vertices[pivot]
    out = {}
    for j in range(g.n):
        if j == pivot:
            continue
        d2 = p.dist2(g.vertices[j])
        if not d2.is_rational():
            continue
        v = d2.c[0]
        if Fraction(1, 4) <= v <= MAXD2:
            out.setdefault(v, []).append(j)
    return out


def sweep(g, label, npivots):
    t = time.time()
    colourable, _ = is_k_colorable(g, K, timeout=600)
    print(f"{label}: {g}  {K}-colourable={colourable}  [{time.time()-t:.0f}s]", flush=True)
    if colourable is False:
        print(f"*** {label} IS ALREADY NOT {K}-COLOURABLE ***", flush=True)
        save_certificate(g, os.path.join(OUT, f"found_k{K}.json"), K,
                         f"chi(R^2) >= {K+1}: no proper {K}-colouring")
        return True
    if colourable is None:
        return False
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:npivots]
    q = 0
    for pivot in order:
        grp = groups_for(g, pivot)
        if not grp:
            continue
        allt = sorted({j for js in grp.values() for j in js})
        test = SeparationTest(g, K, pivot, allt)
        try:
            for val in sorted(grp):
                q += 1
                sep, core = test.run(subset=grp[val])
                if sep:
                    continue
                print(f"  {label} pivot {pivot} d2={val}: FORCED core {len(core)}", flush=True)
                if len(core) != 1:
                    continue
                try:
                    spun, fld = spindle_union_auto(g, pivot, core[0])
                except (ValueError, AssertionError) as e:
                    print(f"    spindle failed: {e}", flush=True)
                    continue
                s2, _ = is_k_colorable(spun, K, timeout=7200)
                print(f"    spindled {spun} over {fld}: {K}-colourable={s2}", flush=True)
                if s2 is False:
                    save_certificate(spun, os.path.join(OUT, f"found_k{K}.json"), K,
                                     f"chi(R^2) >= {K+1}: no proper {K}-colouring")
                    print(f"\n*** NON-{K}-COLOURABLE GRAPH ***", flush=True)
                    return True
        finally:
            test.close()
    print(f"  {label}: {q} queries, nothing forced [{time.time()-t:.0f}s]", flush=True)
    return False


def main():
    g = build_G()
    # tighten by unioning with a rotated copy about the most connected vertex,
    # spindling on a short rational distance that exists there
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    built = 0
    for pivot in order[:40]:
        grp = groups_for(g, pivot)
        for val in sorted(grp):
            if val > 4:
                break
            tgt = grp[val][0]
            try:
                w, fld = spindle_union_auto(g, pivot, tgt)
            except (ValueError, AssertionError):
                continue
            built += 1
            if sweep(w, f"G u rho(G) [pivot {pivot}, d2={val}]", NPIVOTS):
                return
            if built >= 12:
                return
    print("escalation exhausted", flush=True)


if __name__ == "__main__":
    main()
